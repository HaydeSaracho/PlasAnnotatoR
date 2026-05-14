"""
CLI principal de PlasAnnotatoR usando Typer.
"""

import typer
from rich.console import Console
from rich.panel import Panel

from plasannotator import __version__

app = typer.Typer(
    name="plasannotator",
    help="Plasmid prediction, annotation and network analysis tool.",
    add_completion=True,
)

db_app = typer.Typer(help="Gestión de bases de datos.")
app.add_typer(db_app, name="db")

console = Console()


# ---------------------------------------------------------------------------
# Comando raíz
# ---------------------------------------------------------------------------

@app.callback(invoke_without_command=True)
def main(
    ctx: typer.Context,
    version: bool = typer.Option(False, "--version", "-v", help="Muestra la versión."),
) -> None:
    if version:
        console.print(f"PlasAnnotatoR v{__version__}")
        raise typer.Exit()
    if ctx.invoked_subcommand is None:
        console.print(
            Panel.fit(
                f"[bold]PlasAnnotatoR[/bold] v{__version__}\n"
                "Usa [cyan]plasannotator --help[/cyan] para ver los comandos disponibles.",
                border_style="cyan",
            )
        )


# ---------------------------------------------------------------------------
# Subcomandos: db
# ---------------------------------------------------------------------------

@db_app.command("list")
def db_list() -> None:
    """Lista todas las bases de datos y su estado."""
    from plasannotator.db.manager import list_databases
    list_databases()


@db_app.command("index")
def db_index(
    name: str = typer.Argument(..., help="Nombre de la DB (ej: plsdb, refseq)."),
) -> None:
    """Concatena e indexa una base de datos para usarla en predicciones."""
    from plasannotator.db.manager import index_database
    index_database(name)


@db_app.command("index-all")
def db_index_all() -> None:
    """Indexa todas las bases de datos instaladas."""
    from plasannotator.db.manager import index_database, AVAILABLE_DBS
    from plasannotator.config import DB_DIR, AVAILABLE_DBS as DBS

    for name, meta in DBS.items():
        db_path = DB_DIR / meta["dir"]
        if db_path.exists():
            console.print(f"[cyan]Indexando {name}...[/cyan]")
            index_database(name)
        else:
            console.print(f"[yellow]Omitiendo {name} (no instalada)[/yellow]")


# ---------------------------------------------------------------------------
# Subcomandos: predict
# ---------------------------------------------------------------------------

@app.command("predict")
def predict(
    input_fasta: str = typer.Argument(..., help="Archivo FASTA de entrada."),
    db: str = typer.Option("plsdb", "--db", "-d", help="DB a usar (ej: plsdb, refseq)."),
    output: str = typer.Option("results/predictions", "--output", "-o", help="Prefijo de salida."),
    threads: int = typer.Option(4, "--threads", "-t", help="Número de hilos CPU."),
    min_length: int = typer.Option(1000, "--min-length", help="Longitud mínima de contig (pb)."),
    network: bool = typer.Option(False, "--network", "-n", help="Construir red de secuencias."),
) -> None:
    """Predice si los contigs de un FASTA son plasmidios o cromosoma."""
    from pathlib import Path
    from plasannotator.predict.ensemble import run_prediction

    fasta_path = Path(input_fasta)
    if not fasta_path.exists():
        console.print(f"[red]Error:[/red] Archivo no encontrado: {input_fasta}")
        raise typer.Exit(1)

    console.print(f"\n[bold]PlasAnnotatoR[/bold] · predicción")
    console.print(f"  Entrada  : {input_fasta}")
    console.print(f"  DB       : {db}")
    console.print(f"  Hilos    : {threads}")
    console.print(f"  Red      : {'sí' if network else 'no'}\n")

    run_prediction(
        fasta_path=fasta_path,
        db_name=db,
        output_prefix=output,
        threads=threads,
        min_length=min_length,
        build_network=network,
    )


# ---------------------------------------------------------------------------
# Subcomandos: network
# ---------------------------------------------------------------------------

@app.command("network")
def network(
    predictions: str = typer.Argument(..., help="TSV de predicciones previas."),
    db: str = typer.Option("plsdb", "--db", "-d", help="DB de referencia."),
    fasta: str = typer.Option(None, "--fasta", "-f", help="FASTA original de entrada."),
    output: str = typer.Option("results/network", "--output", "-o", help="Prefijo de salida."),
    ani: float = typer.Option(0.95, "--ani", help="Umbral ANI (0-1)."),
    threads: int = typer.Option(4, "--threads", "-t", help="Número de hilos CPU."),
) -> None:
    """Construye una red de secuencias plasmídicas desde un TSV de predicciones."""
    from pathlib import Path
    from plasannotator.network.builder import build_network

    tsv_path = Path(predictions)
    if not tsv_path.exists():
        console.print(f"[red]Error:[/red] Archivo no encontrado: {predictions}")
        raise typer.Exit(1)

    fasta_path = Path(fasta) if fasta else None

    console.print(f"\n[bold]PlasAnnotatoR[/bold] · red de secuencias")
    console.print(f"  Predicciones : {predictions}")
    console.print(f"  DB           : {db}")
    console.print(f"  FASTA        : {fasta or 'auto'}")
    console.print(f"  Umbral ANI   : {ani}\n")

    build_network(
        tsv_path=tsv_path,
        db_name=db,
        output_prefix=output,
        ani_threshold=ani,
        threads=threads,
        original_fasta=fasta_path,
    )


# ---------------------------------------------------------------------------
# Subcomandos: web
# ---------------------------------------------------------------------------

@app.command("web")
def web(
    host: str = typer.Option("127.0.0.1", "--host", help="Host del servidor."),
    port: int = typer.Option(5000, "--port", "-p", help="Puerto del servidor."),
    debug: bool = typer.Option(False, "--debug", help="Modo debug de Flask."),
) -> None:
    """Lanza la interfaz web local de PlasAnnotatoR."""
    from plasannotator.web.app import create_app

    console.print(f"\n[bold]PlasAnnotatoR[/bold] · interfaz web")
    console.print(f"  Abre tu navegador en [cyan]http://{host}:{port}[/cyan]\n")

    flask_app = create_app()
    flask_app.run(host=host, port=port, debug=debug)


if __name__ == "__main__":
    app()