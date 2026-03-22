"""
BiRAGAS CRISPR Complete API Server — FastAPI (DNA + RNA)
"""
import json, logging, os, time, uuid
from typing import Dict, Optional
logger = logging.getLogger("biragas_crispr.api")

try:
    from fastapi import FastAPI, HTTPException, BackgroundTasks
    from fastapi.staticfiles import StaticFiles
    from fastapi.responses import HTMLResponse, JSONResponse
    from fastapi.middleware.cors import CORSMiddleware
    from pydantic import BaseModel
    FASTAPI_OK = True
except ImportError:
    FASTAPI_OK = False

_jobs = {}
_orch = None

def _get_orch():
    global _orch
    if _orch is None:
        from ..pipeline.unified_orchestrator import UnifiedOrchestrator
        _orch = UnifiedOrchestrator({'verbose': True})
    return _orch

def create_app():
    if not FASTAPI_OK:
        raise ImportError("pip install fastapi uvicorn")
    app = FastAPI(title="BiRAGAS CRISPR Complete (DNA + RNA)", version="3.0.0")
    app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

    frontend = os.path.join(os.path.dirname(os.path.dirname(__file__)), "frontend")
    if os.path.isdir(os.path.join(frontend, "static")):
        app.mount("/static", StaticFiles(directory=os.path.join(frontend, "static")), name="static")

    class GuideReq(BaseModel):
        gene_or_sequence: str; n_guides: int = 4; nuclease: str = "NGG"; target_type: str = "auto"
    class AnalyzeReq(BaseModel):
        crispr_dir: str = ""; disease_name: str = "Disease"; output_dir: str = "./output"
        max_knockout_genes: int = 0; run_dna: bool = True; run_rna: bool = True
    class ComboReq(BaseModel):
        gene_a: str; gene_b: str; gene_c: str = ""
    class BaseEditReq(BaseModel):
        gene: str; target_rna: str; edit_type: str = "A-to-I"
    class ModulationReq(BaseModel):
        gene: str; mod_type: str = "CRISPRi"; n_guides: int = 4

    @app.get("/", response_class=HTMLResponse)
    async def index():
        p = os.path.join(frontend, "templates", "index.html")
        if os.path.exists(p):
            with open(p) as f: return f.read()
        # Fallback: serve the standalone app
        standalone = os.path.join(os.path.dirname(os.path.dirname(__file__)), "BiRAGAS_CRISPR_Complete_App.html")
        if os.path.exists(standalone):
            with open(standalone) as f: return f.read()
        return "<h1>BiRAGAS CRISPR Complete</h1>"

    @app.get("/api/status")
    async def status():
        return {"status": "running", "version": "3.0.0", "platform": "DNA + RNA"}

    @app.get("/api/capabilities")
    async def capabilities():
        return _get_orch().get_capabilities()

    @app.post("/api/design-guides")
    async def design_guides(req: GuideReq):
        o = _get_orch(); o._init_engines()
        guides = o._editing.design_guides(req.gene_or_sequence, req.n_guides, req.nuclease, req.target_type)
        return {"gene": req.gene_or_sequence, "nuclease": req.nuclease,
                "target_type": guides[0].target_type if guides else req.target_type,
                "guides": [g.to_dict() for g in guides], "n_returned": len(guides)}

    @app.post("/api/knockout-strategy")
    async def knockout(req: GuideReq):
        o = _get_orch(); o._init_engines()
        s = o._editing.design_knockout_strategy(req.gene_or_sequence, req.n_guides, req.nuclease, req.target_type)
        return {"gene": s.gene, "target_type": s.target_type, "n_configs": s.n_configs,
                "expected_efficiency": s.expected_efficiency,
                "guides": [g.to_dict() for g in s.guides], "configs": s.configs}

    @app.post("/api/base-edit")
    async def base_edit(req: BaseEditReq):
        o = _get_orch(); o._init_engines()
        if req.edit_type == "A-to-I":
            results = o._rna_base_edit.find_best_edit_sites(req.target_rna, "A-to-I", req.gene, 5)
        else:
            results = o._rna_base_edit.find_best_edit_sites(req.target_rna, "C-to-U", req.gene, 5)
        return {"gene": req.gene, "edit_type": req.edit_type,
                "sites": [r.to_dict() for r in results]}

    @app.post("/api/modulation")
    async def modulation(req: ModulationReq):
        o = _get_orch(); o._init_engines()
        if req.mod_type == "CRISPRi":
            guides = o._transcriptome.design_crispri_guides(req.gene, n_guides=req.n_guides)
        else:
            guides = o._transcriptome.design_crispra_guides(req.gene, n_guides=req.n_guides)
        return {"gene": req.gene, "type": req.mod_type,
                "guides": [g.to_dict() for g in guides]}

    @app.post("/api/analyze")
    async def analyze(req: AnalyzeReq, bg: BackgroundTasks):
        jid = str(uuid.uuid4())[:8]
        _jobs[jid] = {'status': 'running', 'started': time.time(), 'result': None}
        def _run():
            try:
                r = _get_orch().run(crispr_dir=req.crispr_dir, disease_name=req.disease_name,
                                     output_dir=req.output_dir, run_dna=req.run_dna, run_rna=req.run_rna)
                _jobs[jid]['result'] = r; _jobs[jid]['status'] = 'completed'
            except Exception as e:
                _jobs[jid]['status'] = 'failed'; _jobs[jid]['error'] = str(e)
        bg.add_task(_run)
        return {"job_id": jid, "status": "running"}

    @app.get("/api/results/{jid}")
    async def results(jid: str):
        if jid not in _jobs: raise HTTPException(404)
        return _jobs[jid]

    @app.get("/api/jobs")
    async def jobs():
        return {k: {'status': v['status']} for k, v in _jobs.items()}

    return app

def run_server(host="0.0.0.0", port=8000):
    if FASTAPI_OK:
        import uvicorn
        print(f"\n{'='*60}\n  BiRAGAS CRISPR Complete v3.0 (DNA + RNA)\n  http://localhost:{port}\n  API docs: http://localhost:{port}/docs\n{'='*60}\n")
        uvicorn.run(create_app(), host=host, port=port)
    else:
        print("Install FastAPI: pip install fastapi uvicorn")
