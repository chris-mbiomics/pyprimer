"""Command-line entry point: ``pyprimer serve`` and quick primer analysis."""

from __future__ import annotations

import argparse
import sys
import webbrowser
from collections.abc import Sequence

from . import __version__
from .report import analyse_pair


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pyprimer",
        description="Primer design for PCR, bisulphite PCR, real-time PCR and sequencing.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    serve = commands.add_parser("serve", help="start the web interface")
    serve.add_argument(
        "--host", default="127.0.0.1", help="interface to bind (default: %(default)s)"
    )
    serve.add_argument("--port", type=int, default=8000, help="port (default: %(default)s)")
    serve.add_argument("--open", action="store_true", help="open the interface in a browser")
    serve.add_argument("--debug", action="store_true", help="run Flask's debug server")

    tm = commands.add_parser("analyse", help="print Tm and primer-dimers for one or two primers")
    tm.add_argument("forward")
    tm.add_argument("reverse", nargs="?")

    args = parser.parse_args(argv)
    if args.command == "serve":
        return _serve(args.host, args.port, args.open, args.debug)
    return _analyse(args.forward, args.reverse)


def _serve(host: str, port: int, open_browser: bool, debug: bool) -> int:
    from .web import create_app

    app = create_app()
    url = f"http://{'localhost' if host in ('127.0.0.1', '0.0.0.0') else host}:{port}/"
    print(f"PyPrimer {__version__} running at {url} (press Ctrl+C to stop)")
    if open_browser:
        webbrowser.open(url)
    if debug:
        app.run(host=host, port=port, debug=True)
        return 0
    try:
        from waitress import serve
    except ImportError:
        app.run(host=host, port=port, threaded=True)
    else:
        serve(app, host=host, port=port)
    return 0


def _analyse(forward: str, reverse: str | None) -> int:
    analysis = analyse_pair(forward, reverse)
    for error in analysis.errors:
        print(f"Error: {error}", file=sys.stderr)
    for name, primer in (("Forward", analysis.forward), ("Reverse", analysis.reverse)):
        if primer:
            print(
                f"{name}: {primer.seq}\n  Tm {primer.tm:.2f} °C, {primer.length} bases, "
                f"{primer.gc:.0f}% GC, dH {primer.dh:.2f} kcal/mol, dS {primer.ds:.2f} eu, "
                f"dG {primer.dg:.2f} kcal/mol"
            )
    print()
    print(analysis.dimer_text())
    return 1 if analysis.errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
