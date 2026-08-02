from __future__ import annotations

import shutil
import subprocess
import sys
import uuid
import zipfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from fastapi import FastAPI, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from tn5cope.input_parsers import parse_sequences_from_text, write_query_csv
from tn5cope.pipeline import TN5_STRUCTURE_PROFILES


PACKAGE_DIR = Path(__file__).resolve().parent
TEMPLATES_DIR = PACKAGE_DIR / "templates"
STATIC_DIR = PACKAGE_DIR / "static"
RUNS_DIR = Path.cwd() / "runs"

RUNS_DIR.mkdir(exist_ok=True)

app = FastAPI(title="Tn5cope")
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))


def safe_filename(name: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name.strip())
    return cleaned or "upload.dat"


def safe_path_component(name: str) -> bool:
    """Accept one ordinary filename component, never dot segments or paths."""

    return bool(
        name
        and name not in {".", ".."}
        and not name.startswith(".")
        and safe_filename(name) == name
    )


def create_job_dir(job_name: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    slug = safe_filename(job_name)[:40] if job_name else "job"
    job_id = f"{stamp}_{slug}_{uuid.uuid4().hex[:8]}"
    job_dir = RUNS_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=False)
    return job_dir


async def save_upload(upload: UploadFile, destination: Path) -> Path:
    with destination.open("wb") as handle:
        while True:
            chunk = await upload.read(1024 * 1024)
            if not chunk:
                break
            handle.write(chunk)
    await upload.close()
    return destination


def bundle_outputs(output_dir: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(output_dir.iterdir()):
            if path.is_file():
                zf.write(path, arcname=path.name)


@app.get("/", response_class=HTMLResponse)
async def index(request: Request) -> HTMLResponse:
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "result": None,
            "error": None,
        },
    )


@app.get("/health")
async def health() -> dict:
    return {"status": "ok"}


@app.post("/analyze", response_class=HTMLResponse)
async def analyze(
    request: Request,
    genome_fasta: UploadFile = File(...),
    gtf_file: UploadFile = File(...),
    operon_file: Optional[UploadFile] = File(None),
    query_file: Optional[UploadFile] = File(None),
    sequences_text: str = Form(""),
    job_name: str = Form("tail_pcr_run"),
    tn5_structure_profile: str = Form("unknown"),
) -> HTMLResponse:
    if tn5_structure_profile not in TN5_STRUCTURE_PROFILES:
        raise HTTPException(status_code=400, detail="Unsupported Tn5 structure profile.")
    sequences: List[str] = []
    sequences.extend(parse_sequences_from_text(sequences_text, "pasted.txt"))

    query_file_name = ""
    if query_file and query_file.filename:
        query_file_name = safe_filename(query_file.filename)
        raw_query_bytes = await query_file.read()
        await query_file.close()
        query_text = raw_query_bytes.decode("utf-8-sig", errors="ignore")
        sequences.extend(parse_sequences_from_text(query_text, query_file_name))

    deduped_sequences = []
    seen = set()
    for sequence in sequences:
        if not sequence or sequence in seen:
            continue
        seen.add(sequence)
        deduped_sequences.append(sequence)

    if not deduped_sequences:
        return templates.TemplateResponse(
            "index.html",
            {
                "request": request,
                "result": None,
                "error": "No valid query sequence was detected. Paste a sequence or upload a CSV/FASTA/TXT query file.",
            },
            status_code=400,
        )

    job_dir = create_job_dir(job_name)
    inputs_dir = job_dir / "inputs"
    outputs_dir = job_dir / "outputs"
    inputs_dir.mkdir()
    outputs_dir.mkdir()

    genome_name = safe_filename(genome_fasta.filename or "reference_genome.fna")
    gtf_name = safe_filename(gtf_file.filename or "reference_annotation.gtf")
    genome_path = await save_upload(genome_fasta, inputs_dir / genome_name)
    gtf_path = await save_upload(gtf_file, inputs_dir / gtf_name)
    operon_name = ""
    operon_path: Optional[Path] = None
    if operon_file and operon_file.filename:
        operon_name = safe_filename(operon_file.filename)
        operon_path = await save_upload(operon_file, inputs_dir / operon_name)
    query_csv_path = inputs_dir / "query.csv"
    write_query_csv(deduped_sequences, query_csv_path)

    output_prefix = "tn5_results"
    command = [
        sys.executable,
        "-m",
        "tn5cope.pipeline",
        "--query-csv",
        str(query_csv_path),
        "--genome-fasta",
        str(genome_path),
        "--gtf",
        str(gtf_path),
        "--output-dir",
        str(outputs_dir),
        "--output-prefix",
        output_prefix,
        "--tn5-structure-profile",
        tn5_structure_profile,
    ]
    if operon_path is not None:
        command.extend(["--operon-tsv", str(operon_path)])
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        shutil.rmtree(job_dir, ignore_errors=True)
        return templates.TemplateResponse(
            "index.html",
            {
                "request": request,
                "result": None,
                "error": f"Analysis failed. stderr: {completed.stderr.strip() or 'unknown error'}",
            },
            status_code=500,
        )

    zip_path = job_dir / f"{output_prefix}_bundle.zip"
    bundle_outputs(outputs_dir, zip_path)

    result = {
        "job_id": job_dir.name,
        "job_name": job_name,
        "sequence_count": len(deduped_sequences),
        "genome_name": genome_name,
        "gtf_name": gtf_name,
        "operon_name": operon_name or "not provided",
        "tn5_structure_profile": tn5_structure_profile,
        "files": [
            path.name
            for path in sorted(outputs_dir.iterdir())
            if path.is_file()
        ],
        "bundle": zip_path.name,
    }
    return templates.TemplateResponse(
        "index.html",
        {
            "request": request,
            "result": result,
            "error": None,
        },
    )


@app.get("/download/{job_id}/{filename}")
async def download(job_id: str, filename: str) -> FileResponse:
    if not safe_path_component(job_id) or not safe_path_component(filename):
        raise HTTPException(status_code=404, detail="File not found.")
    runs_root = RUNS_DIR.resolve()
    job_dir = (runs_root / job_id).resolve()
    if job_dir.parent != runs_root:
        raise HTTPException(status_code=404, detail="Job not found.")
    if not job_dir.exists():
        raise HTTPException(status_code=404, detail="Job not found.")
    candidates = [
        job_dir / filename,
        job_dir / "outputs" / filename,
        job_dir / "inputs" / filename,
    ]
    for candidate in candidates:
        path = candidate.resolve()
        if job_dir not in path.parents:
            continue
        if path.exists() and path.is_file():
            return FileResponse(path)
    raise HTTPException(status_code=404, detail="File not found.")


def run() -> None:
    import uvicorn

    uvicorn.run("tn5cope.web:app", host="127.0.0.1", port=8000, reload=False)
