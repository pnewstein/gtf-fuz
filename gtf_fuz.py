from __future__ import annotations

from pathlib import Path
from os import PathLike
import sys


import numpy as np
from numpy.typing import NDArray

from fuzzywuzzy import process
import click
import polars as pl

DATA_DIR = Path().home() / ".local/share/gtf_fuz"
DATA_DIR.mkdir(exist_ok=True)


def proc_gff(path: PathLike) -> NDArray[np.str_]:
    """
    reads all of the genes (not mt or lncRNA or CG) from a gff
    and returns a list of strings
    """
    cache_path = DATA_DIR / Path(path).with_suffix(".npy").name
    try:
        return np.load(cache_path)
    except FileNotFoundError:
        pass
    from gtfparse import read_gtf
    df = read_gtf(path)
    assert isinstance(df, pl.DataFrame)
    df_genes = (
        df.filter(pl.col("feature") == "gene")
        .filter(
            (~pl.col("gene_name").str.contains(r":|^CG\d+$"))
            & (pl.col("gene_name") != "")
        )
        .select("gene_name")
        .to_series()
        .to_numpy()
        .astype(str)
    )
    np.save(cache_path, df_genes)
    return df_genes


def gtf_fuz(path: PathLike, search: str) -> list[str]:
    """
    returns n_to_show number of genes that are close matches to search string
    """
    gene_names = proc_gff(path)
    names_series = pl.Series(gene_names)
    starts_hit = names_series.filter(names_series.str.to_lowercase().str.starts_with(search.lower()))
    if len(starts_hit) == 0:
        scores = process.extract(search, gene_names)
        return [str(s[0]) for s in scores]
    return starts_hit[starts_hit.str.lengths().arg_sort()].to_list()



@click.command
@click.argument("search-string")
@click.option("--gtf")
def main(search_string: str, gtf: click.Path | None = None):
    if gtf is None:
        try:
            path = Path((DATA_DIR / "last_path").read_text())
        except FileNotFoundError:
            click.echo("Must specify a gtf path", file=sys.stderr)
            sys.exit(4)
    else:
        (DATA_DIR / "last_path").write_text(str(Path(gtf).resolve()))
        path = Path(gtf)
    for string in gtf_fuz(path, search_string):
        click.echo(f'"{string}",')

