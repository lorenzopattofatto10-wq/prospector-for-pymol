"""One-time preparation (not shipped): straight, ready-to-pack shapes of lipids
and detergents from the PDB dictionary's ideal geometry.  Bond lengths and
angles are kept; only torsions are changed:
  - every C-C-C-C torsion along a tail is set to trans (180 deg); a cis
    double bond stays cis (the real kink of an oleoyl chain)
  - ester C-O bonds are set to trans (their normal state)
  - the torsions between head and tails are searched (60/180/300 deg) for the
    slimmest shape with both tails parallel and no atoms bumping
Output: heavy atoms only, head atom at the origin, tails along -z, plus the
tightest clash-free staggered lattice for packing them side by side."""
import itertools, json, math
import numpy as np

def read_sdf(path):
    lines = open(path).read().split("\n")
    na, nb = int(lines[3][:3]), int(lines[3][3:6])
    atoms = [(l[31:34].strip(), [float(l[0:10]), float(l[10:20]), float(l[20:30])])
             for l in lines[4:4 + na]]
    bonds = [(int(l[0:3]) - 1, int(l[3:6]) - 1, int(l[6:9])) for l in lines[4 + na:4 + na + nb]]
    keep = [i for i, a in enumerate(atoms) if a[0] != "H"]
    new = {old: n for n, old in enumerate(keep)}
    return ([atoms[i][0] for i in keep], np.array([atoms[i][1] for i in keep]),
            [(new[a], new[b], o) for a, b, o in bonds if a in new and b in new])

def write_mol(path, title, elem, xyz, bonds):
    out = [title, "  prospector straightened conformer", "",
           "%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (len(elem), len(bonds))]
    for e, (x, y, z) in zip(elem, xyz):
        out.append("%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0" % (x, y, z, e))
    for a, b, o in bonds:
        out.append("%3d%3d%3d  0" % (a + 1, b + 1, o))
    out.append("M  END")
    open(path, "w").write("\n".join(out) + "\n$$$$\n")

def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n; w = b2 - np.dot(b2, b1n) * b1n
    return math.degrees(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))

def rotate(xyz, idx, origin, axis, angle):
    axis = axis / np.linalg.norm(axis); a = math.radians(angle)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    R = np.eye(3) + math.sin(a) * K + (1 - math.cos(a)) * K @ K
    xyz[idx] = (xyz[idx] - origin) @ R.T + origin

def prepare(code, head_element, tails):
    elem, xyz, bonds = read_sdf(code + ".sdf")
    n = len(elem)
    nb = {i: [] for i in range(n)}; order = {}
    for a, b, o in bonds:
        nb[a].append(b); nb[b].append(a); order[(a, b)] = order[(b, a)] = o
    def bfs(start, banned=None):
        depth = {start: 0}; prev = {start: None}; queue = [start]
        while queue:
            x = queue.pop(0)
            for y in nb[x]:
                if y not in depth and (banned is None or (x, y) not in banned):
                    depth[y] = depth[x] + 1; prev[y] = x; queue.append(y)
        return depth, prev
    head = [i for i in range(n) if elem[i] == head_element]
    if head_element == "O":       # the hydroxyl / sugar oxygen farthest from the tail
        terminal = [i for i in range(n) if elem[i] == "C" and len(nb[i]) == 1]
        dmax = {i: max(bfs(i)[0].values()) for i in head}
        head = sorted(head, key=lambda i: -max(bfs(i)[0][t] for t in terminal))
    head = head[0]
    depth, prev = bfs(head)
    terminal = [i for i in range(n) if elem[i] == "C" and len(nb[i]) == 1]
    ends = sorted(terminal, key=lambda i: -depth[i])[:tails]
    routes = []
    for e in ends:
        r, x = [], e
        while x is not None: r.append(x); x = prev[x]
        routes.append(r[::-1])
    def ring(b, c):
        d, _ = bfs(b, banned={(b, c), (c, b)})
        return c in d
    def side(b, c):            # atoms on c's side of bond b-c
        d, _ = bfs(c, banned={(c, b), (b, c)})
        return sorted(d)
    def set_t(xyz, q, value):
        a, b, c, d = q
        rotate(xyz, side(b, c), xyz[b], xyz[c] - xyz[b], value - dihedral(*xyz[[a, b, c, d]]))
    tail_q, ester_q, linker_q = [], [], []
    for r in routes:
        for k in range(len(r) - 3):
            a, b, c, d = r[k:k + 4]
            if order[(b, c)] != 1 or ring(b, c):
                continue
            q = (a, b, c, d)
            if all(elem[x] == "C" for x in r[k:]):
                if q not in tail_q: tail_q.append(q)
            elif elem[b] == "O" and elem[c] == "C" and any(elem[x] == "O" and order[(c, x)] == 2
                                                           for x in nb[c]):
                if q not in ester_q: ester_q.append(q)
            elif all(set((b, c)) != set((x[1], x[2])) for x in linker_q):
                linker_q.append(q)
    for q in tail_q + ester_q:
        set_t(xyz, q, 180.0)
    # pairs more than 3 bonds apart, for the bump check
    sep = np.array([[bfs(i)[0].get(j, 99) for j in range(n)] for i in range(n)])
    far = sep > 4                      # more than 4 bonds apart: at least 3.0 A
    near = sep == 4                    # 1-5 neighbours sit closer by nature: 2.5 A
    best = None
    linker_q = linker_q[:9]
    for combo in itertools.product((60.0, 180.0, 300.0), repeat=len(linker_q)):
        trial = xyz.copy()
        for q, v in zip(linker_q, combo):
            set_t(trial, q, v)
        d = np.linalg.norm(trial[:, None, :] - trial[None, :, :], axis=2)
        if (d[far] < 3.0).any() or (d[near] < 2.5).any():
            continue
        t = trial[ends].mean(axis=0); axis = t - trial[head]; L = np.linalg.norm(axis); u = axis / L
        v = trial - trial[head]; along = v @ u
        perp = np.sqrt(np.maximum(0, (v * v).sum(1) - along ** 2))
        spread = np.linalg.norm(trial[ends[0]] - trial[ends[-1]]) if tails > 1 else 0
        score = perp.max() + 0.6 * spread - 0.3 * L
        if best is None or score < best[0]:
            best = (score, trial.copy(), L)
    xyz = best[1]
    # head at the origin, tails down -z, widest direction along x
    xyz = xyz - xyz[head]
    u = xyz[ends].mean(axis=0); u = u / np.linalg.norm(u); w = np.array([0, 0, -1.0])
    ax = np.cross(u, w)
    if np.linalg.norm(ax) > 1e-6:
        rotate(xyz, list(range(n)), np.zeros(3), ax, math.degrees(math.atan2(np.linalg.norm(ax), u @ w)))
    best_angle = min(range(0, 180, 3), key=lambda g: np.ptp(
        (xyz[:, :2] @ np.array([[math.cos(math.radians(g)), -math.sin(math.radians(g))],
                                [math.sin(math.radians(g)), math.cos(math.radians(g))]]).T)[:, 1]))
    rotate(xyz, list(range(n)), np.zeros(3), np.array([0, 0, 1.0]), -best_angle)
    if np.ptp(xyz[:, 1]) > np.ptp(xyz[:, 0]):
        rotate(xyz, list(range(n)), np.zeros(3), np.array([0, 0, 1.0]), 90)
    # tightest clash-free staggered lattice (3.0 A between molecules)
    reach = float(max(np.ptp(xyz[:, 0]), np.ptp(xyz[:, 1]))) + 3.0
    def clash_free(sx, sy):
        shifts = []                      # every lattice neighbour close enough to touch
        for j in range(0, int(reach // sy) + 2):
            for i in range(-int(reach // sx) - 2, int(reach // sx) + 3):
                dx, dy = i * sx + (sx / 2 if j % 2 else 0.0), j * sy
                if (j > 0 or dx > 0) and math.hypot(dx, dy) < reach:
                    shifts.append((dx, dy))
        for dx, dy in shifts:
            for dz in (0.0, 1.0, -1.0):
                other = xyz + np.array([dx, dy, dz])
                if (np.linalg.norm(xyz[:, None] - other[None], axis=2) < 3.0).any():
                    return False
        return True
    lattice = None
    for sx in np.arange(4.0, 25.0, 0.25):
        for sy in np.arange(4.0, 25.0, 0.25):
            if lattice and sx * sy >= lattice[0]: continue
            if clash_free(sx, sy): lattice = (sx * sy, sx, sy)
    write_mol(code + ".mol", code, elem, xyz, bonds)
    info = {"length": round(float(-xyz[:, 2].min()), 2), "lattice": [lattice[1], lattice[2]],
            "area": round(lattice[0], 1), "atoms": n,
            "width": [round(float(np.ptp(xyz[:, 0])), 1), round(float(np.ptp(xyz[:, 1])), 1)]}
    print("%-4s %d tail + %d ester + %d linker torsions | length %.1f A | cross-section %s A | "
          "lattice %.2f x %.2f A = %.0f A^2 per molecule"
          % (code, len(tail_q), len(ester_q), len(linker_q), info["length"], info["width"],
             lattice[1], lattice[2], lattice[0]))
    return info

table = {}
for code, head, tails in (("POV", "P", 2), ("PCW", "P", 2), ("PEE", "P", 2), ("CLR", "O", 1),
                          ("SDS", "S", 1), ("LDA", "N", 1)):
    table[code] = prepare(code, head, tails)
json.dump(table, open("shapes.json", "w"), indent=1)
