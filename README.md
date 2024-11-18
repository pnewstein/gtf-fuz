# gtf-fuz

A fuzzy searcher for gene names obtained from a gtf-file

It first returns all genes whose name starts with the search (case insensitive)

If none exist, it uses a fuzzy finder to find

## usage
```bash
gtf-fuz --gtf path.gtf hb
```

After your first run, it caches the gene names and --gtf is no longer necessary

```bash
gtf-fuz --gtf path.gtf hb
```

## installation
pip install git:https://www.github.com/pnewstein/gtf-fuz
