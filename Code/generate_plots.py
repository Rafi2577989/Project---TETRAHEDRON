from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
WALL_CSV = BASE / "wall_experiment.csv"
SPHERE_CSV = BASE / "sphere_experiment.csv"
OUT = BASE / "plots"
OUT.mkdir(exist_ok=True)

for old_pdf in OUT.glob("*.pdf"):
    old_pdf.unlink()

REST_LENGTH = 1.0
SPHERE_RADIUS = 0.75
SPRING_PAIRS = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]

wall = pd.read_csv(WALL_CSV)
sphere = pd.read_csv(SPHERE_CSV)

def save_current(name):
    plt.tight_layout()
    plt.savefig(OUT / f"{name}.png", dpi=240, bbox_inches="tight")
    plt.close()

def spring_deformations(df):
    result = {}
    for i, j in SPRING_PAIRS:
        dx = df[f"p{j}x"] - df[f"p{i}x"]
        dy = df[f"p{j}y"] - df[f"p{i}y"]
        dz = df[f"p{j}z"] - df[f"p{i}z"]
        L = np.sqrt(dx*dx + dy*dy + dz*dz)
        result[f"P{i}-P{j}"] = L - REST_LENGTH
    return pd.DataFrame(result)

# 1. Wall energy transfer
plt.figure(figsize=(10, 6))
plt.plot(wall["time"], wall["com_kinetic"], label="COM kinetic energy")
plt.plot(wall["time"], wall["vibration_total"], label="Vibration energy")
plt.plot(wall["time"], wall["total"], label="Total energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("1. Wall collision: energy transfer")
plt.legend()
plt.grid(True, alpha=0.25)
save_current("01_wall_energy_transfer")

# 2. Sphere energy transfer
plt.figure(figsize=(10, 6))
plt.plot(sphere["time"], sphere["com_kinetic"], label="COM kinetic energy")
plt.plot(sphere["time"], sphere["rotational"], label="Rotational energy")
plt.plot(sphere["time"], sphere["vibration_total"], label="Vibration energy")
plt.plot(sphere["time"], sphere["total"], label="Total energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("2. Sphere collision: energy transfer")
plt.legend()
plt.grid(True, alpha=0.25)
save_current("02_sphere_energy_transfer")

# 3. Wall kinetic vs spring energy
plt.figure(figsize=(10, 6))
plt.plot(wall["time"], wall["kinetic"], label="Total kinetic energy")
plt.plot(wall["time"], wall["spring"], label="Spring potential energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("3. Wall collision: kinetic and spring energy")
plt.legend()
plt.grid(True, alpha=0.25)
save_current("03_wall_kinetic_vs_spring")

# 4. Sphere kinetic vs spring energy
plt.figure(figsize=(10, 6))
plt.plot(sphere["time"], sphere["kinetic"], label="Total kinetic energy")
plt.plot(sphere["time"], sphere["spring"], label="Spring potential energy")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.title("4. Sphere collision: kinetic and spring energy")
plt.legend()
plt.grid(True, alpha=0.25)
save_current("04_sphere_kinetic_vs_spring")

# 5. Normalized total energy comparison
plt.figure(figsize=(10, 6))
plt.plot(wall["time"], wall["total"] / wall["total"].iloc[0], label="Wall: E/E0")
plt.plot(sphere["time"], sphere["total"] / sphere["total"].iloc[0], label="Sphere: E/E0")
plt.xlabel("Time")
plt.ylabel("Normalized total energy")
plt.title("5. Total-energy conservation: wall vs sphere")
plt.legend()
plt.grid(True, alpha=0.25)
save_current("05_total_energy_conservation_comparison")

# 6. Relative energy error comparison
plt.figure(figsize=(10, 6))
plt.plot(wall["time"], wall["relative_energy_error"], label="Wall")
plt.plot(sphere["time"], sphere["relative_energy_error"], label="Sphere")
plt.xlabel("Time")
plt.ylabel("Relative energy error")
plt.title("6. Numerical energy error: wall vs sphere")
plt.legend()
plt.grid(True, alpha=0.25)
save_current("06_relative_energy_error_comparison")

# 7. Wall COM 3D trajectory
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")
ax.plot(wall["com_x"], wall["com_y"], wall["com_z"], label="COM path")
ax.scatter([wall["com_x"].iloc[0]], [wall["com_y"].iloc[0]], [wall["com_z"].iloc[0]], label="Start")
ax.scatter([wall["com_x"].iloc[-1]], [wall["com_y"].iloc[-1]], [wall["com_z"].iloc[-1]], label="End")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("7. Wall experiment: centre-of-mass trajectory")
ax.legend()
save_current("07_wall_com_trajectory_3d")

# 8. Sphere COM 3D trajectory
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")
ax.plot(sphere["com_x"], sphere["com_y"], sphere["com_z"], label="COM path")
ax.scatter([sphere["com_x"].iloc[0]], [sphere["com_y"].iloc[0]], [sphere["com_z"].iloc[0]], label="Start")
ax.scatter([sphere["com_x"].iloc[-1]], [sphere["com_y"].iloc[-1]], [sphere["com_z"].iloc[-1]], label="End")

u = np.linspace(0, 2*np.pi, 30)
v = np.linspace(0, np.pi, 18)
xs = SPHERE_RADIUS * np.outer(np.cos(u), np.sin(v))
ys = SPHERE_RADIUS * np.outer(np.sin(u), np.sin(v))
zs = SPHERE_RADIUS * np.outer(np.ones_like(u), np.cos(v))
ax.plot_wireframe(xs, ys, zs, rstride=3, cstride=3, linewidth=0.5)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("8. Sphere experiment: centre-of-mass trajectory")
ax.legend()
save_current("08_sphere_com_trajectory_3d")

# 9. Wall spring deformation
wall_def = spring_deformations(wall)
plt.figure(figsize=(10, 6))
for col in wall_def.columns:
    plt.plot(wall["time"], wall_def[col], label=col)
plt.axhline(0, linewidth=1)
plt.xlabel("Time")
plt.ylabel("Spring deformation L - a")
plt.title("9. Wall collision: deformation of all six springs")
plt.legend(ncol=2)
plt.grid(True, alpha=0.25)
save_current("09_wall_spring_deformation")

# 10. Sphere spring deformation
sphere_def = spring_deformations(sphere)
plt.figure(figsize=(10, 6))
for col in sphere_def.columns:
    plt.plot(sphere["time"], sphere_def[col], label=col)
plt.axhline(0, linewidth=1)
plt.xlabel("Time")
plt.ylabel("Spring deformation L - a")
plt.title("10. Sphere collision: deformation of all six springs")
plt.legend(ncol=2)
plt.grid(True, alpha=0.25)
save_current("10_sphere_spring_deformation")

print(f"Done. PNG plots saved in: {OUT}")