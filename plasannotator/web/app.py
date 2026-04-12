"""
Interfaz web de PlasAnnotatoR usando Flask.
"""

from pathlib import Path

from flask import Flask, jsonify, render_template, request

from plasannotator.config import AVAILABLE_DBS, DB_DIR, INDEX_DIR, RESULTS_DIR, logger


def create_app() -> Flask:
    """Factory de la aplicación Flask."""
    app = Flask(__name__, template_folder="templates", static_folder="static")

    # -----------------------------------------------------------------------
    # Rutas principales
    # -----------------------------------------------------------------------

    @app.route("/")
    def index():
        return render_template("index.html")

    # -----------------------------------------------------------------------
    # API: bases de datos
    # -----------------------------------------------------------------------

    @app.route("/api/db/list")
    def api_db_list():
        dbs = []
        for name, meta in AVAILABLE_DBS.items():
            db_path = DB_DIR / meta["dir"]
            idx_path = INDEX_DIR / f"{name}.mmi"
            dbs.append({
                "name": name,
                "description": meta["description"],
                "sequences": meta["sequences"],
                "installed": db_path.exists(),
                "indexed": idx_path.exists(),
            })
        return jsonify(dbs)

    @app.route("/api/db/index/<name>", methods=["POST"])
    def api_db_index(name: str):
        from plasannotator.db.manager import index_database
        try:
            index_database(name)
            return jsonify({"status": "ok", "message": f"'{name}' indexada correctamente."})
        except Exception as e:
            return jsonify({"status": "error", "message": str(e)}), 400

    # -----------------------------------------------------------------------
    # API: predicción
    # -----------------------------------------------------------------------

    @app.route("/api/predict", methods=["POST"])
    def api_predict():
        from plasannotator.predict.ensemble import run_prediction

        if "fasta" not in request.files:
            return jsonify({"status": "error", "message": "No se recibió archivo FASTA."}), 400

        fasta_file = request.files["fasta"]
        db_name = request.form.get("db", "plsdb")
        threads = int(request.form.get("threads", 4))
        min_length = int(request.form.get("min_length", 1000))
        build_network = request.form.get("network", "false").lower() == "true"

        # Guardar FASTA temporal
        upload_dir = RESULTS_DIR / "uploads"
        upload_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = upload_dir / fasta_file.filename
        fasta_file.save(str(fasta_path))

        output_prefix = str(RESULTS_DIR / fasta_path.stem / "predictions")

        try:
            run_prediction(
                fasta_path=fasta_path,
                db_name=db_name,
                output_prefix=output_prefix,
                threads=threads,
                min_length=min_length,
                build_network=build_network,
            )
            tsv_path = Path(f"{output_prefix}.tsv")
            results = _read_tsv(tsv_path)
            return jsonify({"status": "ok", "results": results})
        except Exception as e:
            logger.error(str(e))
            return jsonify({"status": "error", "message": str(e)}), 500

    # -----------------------------------------------------------------------
    # API: red
    # -----------------------------------------------------------------------

    @app.route("/api/network", methods=["POST"])
    def api_network():
        from plasannotator.network.builder import build_network

        data = request.get_json()
        tsv_path = Path(data.get("tsv_path", ""))
        db_name = data.get("db", "plsdb")
        ani_threshold = float(data.get("ani", 0.95))
        threads = int(data.get("threads", 4))

        if not tsv_path.exists():
            return jsonify({"status": "error", "message": f"TSV no encontrado: {tsv_path}"}), 400

        output_prefix = str(tsv_path.parent / "network")

        try:
            build_network(
                tsv_path=tsv_path,
                db_name=db_name,
                output_prefix=output_prefix,
                ani_threshold=ani_threshold,
                threads=threads,
            )
            json_path = Path(f"{output_prefix}_network.json")
            return jsonify({
                "status": "ok",
                "network_json": str(json_path),
            })
        except Exception as e:
            logger.error(str(e))
            return jsonify({"status": "error", "message": str(e)}), 500

    # -----------------------------------------------------------------------
    # API: resultados
    # -----------------------------------------------------------------------

    @app.route("/api/results")
    def api_results():
        results = []
        for tsv in RESULTS_DIR.rglob("predictions.tsv"):
            results.append({
                "name": tsv.parent.name,
                "path": str(tsv),
            })
        return jsonify(results)

    return app


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def _read_tsv(tsv_path: Path) -> list[dict]:
    """Lee un TSV de predicciones y retorna lista de dicts."""
    import csv
    rows = []
    with open(tsv_path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows