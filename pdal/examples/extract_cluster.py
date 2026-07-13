"""Extract a decimated, reprojected point set from ALL .laz/.las tiles in a dir
(a 3x3 OS-grid cluster drops in unchanged). One PDAL pipeline: reader per tile
-> merge -> decimate -> reproject 27700->4326 -> CSV (X=lon,Y=lat,Z,Classification).

Usage: extract_cluster.py <tile_dir> <out.csv> [step]
Needs the `pdal` CLI on PATH.
"""
import sys, os, json, glob, subprocess

tile_dir, out_csv = sys.argv[1], sys.argv[2]
step = int(sys.argv[3]) if len(sys.argv) > 3 else 200
tiles = sorted(glob.glob(os.path.join(tile_dir, '*.laz')) +
               glob.glob(os.path.join(tile_dir, '*.las')))
if not tiles:
    sys.exit(f'no .laz/.las in {tile_dir}')

stages = [{"type": "readers.las", "filename": t, "tag": f"t{i}"}
          for i, t in enumerate(tiles)]
stages.append({"type": "filters.merge", "inputs": [f"t{i}" for i in range(len(tiles))]})
stages.append({"type": "filters.decimation", "step": step})
stages.append({"type": "filters.reprojection",
               "in_srs": "EPSG:27700", "out_srs": "EPSG:4326"})
stages.append({"type": "writers.text", "filename": out_csv,
               "order": "X,Y,Z,Classification",
               "keep_unspecified": "false", "precision": 9})

pipe = os.path.join(os.path.dirname(os.path.abspath(out_csv)), "cluster_pipeline.json")
with open(pipe, "w") as f:
    json.dump(stages, f, indent=2)
print(f'{len(tiles)} tile(s): {", ".join(os.path.basename(t) for t in tiles)}')
subprocess.run(["pdal", "pipeline", pipe], check=True,
               env=dict(os.environ, OMP_NUM_THREADS="8"))
n = sum(1 for _ in open(out_csv)) - 1
print(f'wrote {n} points (1/{step}) -> {out_csv}')
