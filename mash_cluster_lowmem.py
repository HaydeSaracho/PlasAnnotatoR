import sys

def main():
    infile = sys.argv[1]
    threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    outfile = sys.argv[3] if len(sys.argv) > 3 else "representatives.txt"

    print("Pasada 1: leyendo IDs...")
    ids = []
    with open(infile) as f:
        n = int(f.readline().strip())
        for line in f:
            ids.append(line.split('\t')[0])
    print(f"  {len(ids)} secuencias")

    assigned = set()
    representatives = []

    print(f"Pasada 2: clustering greedy (umbral={threshold})...")
    with open(infile) as f:
        f.readline()
        for i, line in enumerate(f):
            parts = line.strip().split('\t')
            seq_id = parts[0]
            dists = list(map(float, parts[1:])) if len(parts) > 1 else []

            if seq_id not in assigned:
                representatives.append(seq_id)
                assigned.add(seq_id)
                for j, d in enumerate(dists):
                    if d <= threshold:
                        assigned.add(ids[j])

            if i % 5000 == 0:
                print(f"  {i}/{len(ids)} · {len(representatives)} representantes")

    with open(outfile, 'w') as f:
        for r in representatives:
            f.write(r + '\n')

    print(f"  {len(representatives)} representantes → {outfile}")

if __name__ == '__main__':
    main()
