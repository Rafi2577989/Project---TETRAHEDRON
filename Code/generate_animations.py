import matplotlib
matplotlib.use('Agg')

from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.gridspec import GridSpec

BASE = Path(__file__).resolve().parent
WALL_CSV = BASE / 'wall_experiment.csv'
SPHERE_CSV = BASE / 'sphere_experiment.csv'
ANIMDIR = BASE / 'animations'
ANIMDIR.mkdir(exist_ok=True)

PAIRS = [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
SPHERE_RADIUS = 0.75
FPS = 12

wall = pd.read_csv(WALL_CSV)
sphere = pd.read_csv(SPHERE_CSV)


def sample_indices(df, max_frames):
    return np.linspace(0, len(df)-1, min(max_frames, len(df)), dtype=int)


def draw_box(ax, bmin=-3.0, bmax=3.0):
    corners = np.array([
        [bmin,bmin,bmin],[bmax,bmin,bmin],[bmax,bmax,bmin],[bmin,bmax,bmin],
        [bmin,bmin,bmax],[bmax,bmin,bmax],[bmax,bmax,bmax],[bmin,bmax,bmax]
    ])
    edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    for a,b in edges:
        ax.plot([corners[a,0], corners[b,0]],
                [corners[a,1], corners[b,1]],
                [corners[a,2], corners[b,2]], linewidth=0.7)


def draw_sphere(ax, radius=SPHERE_RADIUS):
    u = np.linspace(0, 2*np.pi, 16)
    v = np.linspace(0, np.pi, 10)
    xs = radius * np.outer(np.cos(u), np.sin(v))
    ys = radius * np.outer(np.sin(u), np.sin(v))
    zs = radius * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, rstride=2, cstride=2, linewidth=0.35)


def save_animation(anim, path):
    writer = FFMpegWriter(fps=FPS, bitrate=1200)
    anim.save(str(path), writer=writer, dpi=90)


def xyz_from_row(row):
    return np.array([[row[f'p{i}x'], row[f'p{i}y'], row[f'p{i}z']] for i in range(4)])


def make_motion_animation(df, output_path, title, with_sphere=False, max_frames=100):
    idx = sample_indices(df, max_frames)
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-3.2, 3.2)
    ax.set_ylim(-3.2, 3.2)
    ax.set_zlim(-3.2, 3.2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    ax.view_init(elev=24, azim=35)
    draw_box(ax)
    if with_sphere:
        draw_sphere(ax)

    points = ax.scatter([], [], [], s=32)
    spring_lines = [ax.plot([], [], [], linewidth=1.2)[0] for _ in PAIRS]
    com_trail = ax.plot([], [], [], linewidth=1.0, label='COM trail')[0]
    ax.legend(loc='upper left', fontsize=8)
    time_text = ax.text2D(0.02, 0.95, '', transform=ax.transAxes)
    history = []

    def init():
        points._offsets3d = ([], [], [])
        for line in spring_lines + [com_trail]:
            line.set_data([], [])
            line.set_3d_properties([])
        time_text.set_text('')
        return [points, *spring_lines, com_trail, time_text]

    def update(frame_no):
        row = df.iloc[idx[frame_no]]
        xyz = xyz_from_row(row)
        points._offsets3d = (xyz[:,0], xyz[:,1], xyz[:,2])

        for line, (i, j) in zip(spring_lines, PAIRS):
            line.set_data([xyz[i,0], xyz[j,0]], [xyz[i,1], xyz[j,1]])
            line.set_3d_properties([xyz[i,2], xyz[j,2]])

        history.append((row['com_x'], row['com_y'], row['com_z']))
        arr = np.array(history[-45:])
        com_trail.set_data(arr[:,0], arr[:,1])
        com_trail.set_3d_properties(arr[:,2])
        time_text.set_text(f"t = {row['time']:.3f}")
        return [points, *spring_lines, com_trail, time_text]

    anim = FuncAnimation(fig, update, frames=len(idx), init_func=init,
                         interval=1000/FPS, blit=False)
    save_animation(anim, output_path)
    plt.close(fig)


def make_combined_motion_and_energy(df, output_path, title, with_sphere=True, max_frames=80):
    idx = sample_indices(df, max_frames)
    fig = plt.figure(figsize=(10, 5.2))
    gs = GridSpec(1, 2, width_ratios=[1.0, 1.0], figure=fig)
    ax3 = fig.add_subplot(gs[0,0], projection='3d')
    ax2 = fig.add_subplot(gs[0,1])

    ax3.set_xlim(-3.2, 3.2)
    ax3.set_ylim(-3.2, 3.2)
    ax3.set_zlim(-3.2, 3.2)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_zlabel('z')
    ax3.set_title('Motion')
    ax3.view_init(elev=24, azim=35)
    draw_box(ax3)
    if with_sphere:
        draw_sphere(ax3)

    points = ax3.scatter([], [], [], s=32)
    spring_lines = [ax3.plot([], [], [], linewidth=1.2)[0] for _ in PAIRS]
    com_trail = ax3.plot([], [], [], linewidth=1.0)[0]
    time_text = ax3.text2D(0.02, 0.95, '', transform=ax3.transAxes)
    history = []

    ax2.set_xlim(df['time'].min(), df['time'].max())
    ymax = max(df['total'].max(), df['com_kinetic'].max(), df['vibration_total'].max(), df['rotational'].max()) * 1.08
    ax2.set_ylim(0, ymax)
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Energy')
    ax2.set_title('Live energies')
    ax2.grid(True, alpha=0.25)

    line_com = ax2.plot([], [], label='COM kinetic')[0]
    line_rot = ax2.plot([], [], label='Rotation')[0]
    line_vib = ax2.plot([], [], label='Vibration')[0]
    line_tot = ax2.plot([], [], label='Total')[0]
    marker = ax2.axvline(df['time'].iloc[0], linestyle='--')
    ax2.legend(fontsize=8)

    def init():
        points._offsets3d = ([], [], [])
        for line in spring_lines + [com_trail]:
            line.set_data([], [])
            line.set_3d_properties([])
        for line in [line_com, line_rot, line_vib, line_tot]:
            line.set_data([], [])
        marker.set_xdata([df['time'].iloc[0], df['time'].iloc[0]])
        time_text.set_text('')
        return [points, *spring_lines, com_trail, line_com, line_rot, line_vib, line_tot, marker, time_text]

    def update(frame_no):
        ii = idx[frame_no]
        row = df.iloc[ii]
        xyz = xyz_from_row(row)
        points._offsets3d = (xyz[:,0], xyz[:,1], xyz[:,2])

        for line, (i, j) in zip(spring_lines, PAIRS):
            line.set_data([xyz[i,0], xyz[j,0]], [xyz[i,1], xyz[j,1]])
            line.set_3d_properties([xyz[i,2], xyz[j,2]])

        history.append((row['com_x'], row['com_y'], row['com_z']))
        arr = np.array(history[-45:])
        com_trail.set_data(arr[:,0], arr[:,1])
        com_trail.set_3d_properties(arr[:,2])
        time_text.set_text(f"t = {row['time']:.3f}")

        end = ii + 1
        t = df['time'].iloc[:end]
        line_com.set_data(t, df['com_kinetic'].iloc[:end])
        line_rot.set_data(t, df['rotational'].iloc[:end])
        line_vib.set_data(t, df['vibration_total'].iloc[:end])
        line_tot.set_data(t, df['total'].iloc[:end])
        x = df['time'].iloc[ii]
        marker.set_xdata([x, x])
        return [points, *spring_lines, com_trail, line_com, line_rot, line_vib, line_tot, marker, time_text]

    fig.suptitle(title)
    anim = FuncAnimation(fig, update, frames=len(idx), init_func=init,
                         interval=1000/FPS, blit=False)
    save_animation(anim, output_path)
    plt.close(fig)


# Create only 3 videos
make_motion_animation(
    wall,
    ANIMDIR / '01_wall_collision_motion.mp4',
    '1. Tetrahedron wall collision',
    with_sphere=False,
    max_frames=100
)

make_motion_animation(
    sphere,
    ANIMDIR / '02_sphere_collision_motion.mp4',
    '2. Tetrahedron central-sphere collision',
    with_sphere=True,
    max_frames=100
)

make_combined_motion_and_energy(
    sphere,
    ANIMDIR / '03_combined_motion_and_energy_sphere.mp4',
    '6. Motion with live energy readout',
    with_sphere=True,
    max_frames=80
)

print('Done. Created folder:')
print(ANIMDIR)
print('Created files:')
print('01_wall_collision_motion.mp4')
print('02_sphere_collision_motion.mp4')
print('03_combined_motion_and_energy_sphere.mp4')