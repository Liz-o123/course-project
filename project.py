from vpython import *

# Get the temperature input from the user
T = float(input("Enter the temperature (in Kelvin): "))

# Constants and initial setup
k = 1.38e-23  # Boltzmann constant
N = 60  # Number of helium atoms
L = ((24.2e-3 / (6e23)) * N) ** (1 / 3.0) / 2  # Calculates half the side length of a cubic container that can hold 60 helium atoms, based on the molar volume of helium gas at standard temperature and pressure.
m, size = 4e-3 / (6e23), 310e-12  # Mass and size of helium atoms
L_size = L - size # Adjusting the container size to avoid overlap with molecule size.
t, dt = 0, 0.5e-13 # t is the initial time, dt is the time step (0.5E-13 seconds).
vrms = (3 * k * T / m) ** 0.5  # Calculate vrms based on the entered temperature


#Initialization
scene = canvas(width=800, height=800, background=vector(0.1, 0.2, 0.1))

container = box(length=2 * L, height=2 * L, width=2 * L, opacity=0.2, color=color.yellow)

# Initialize the pressure label
pressure_label = label(pos=vector(0, 1.1*L_size, 0), text='Pressure: 0 Pa', height=20, box=False, color=color.white)

atoms = []
momentum = 0  # Initialize momentum here

for i in range(N):
    position = vector(-L_size + 2 * L_size * random(), -L_size + 2 * L_size * random(), -L_size + 2 * L_size * random())

    if i == N - 1:
        atom = sphere(pos=position, radius=size, color=color.yellow, make_trail=True, retain=600)
    else:
        atom = sphere(pos=position, radius=size, color=vector(random(), random(), random()))

    ra, rb = pi * random(), 2 * pi * random()
    atom.m, atom.v = m, vector(vrms * sin(ra) * cos(rb), vrms * sin(ra) * sin(rb), vrms * cos(ra))
    atoms.append(atom)

# Function to handle elastic collisions between atoms
def vcollision(a1, a2):
    v1prime = (a1.v - 2 * a2.m / (a1.m + a2.m) * (a1.pos - a2.pos) * dot(a1.v - a2.v, a1.pos - a2.pos) / mag(a1.pos - a2.pos) ** 2)
    v2prime = (a2.v - 2 * a1.m / (a1.m + a2.m) * (a2.pos - a1.pos) * dot(a2.v - a1.v, a2.pos - a1.pos) / mag(a2.pos - a1.pos) ** 2)
    return v1prime, v2prime

while True:
    t += dt  # Increment the simulation time by the time step

    rate(1000)

    # Calculate new positions for all atoms based on their current velocities
    for i in range(N):
        atoms[i].pos += atoms[i].v * dt

    # Find and handle collisions between pairs of atoms
    for i in range(N - 1):
        for j in range(i + 1, N):
            if mag(atoms[i].pos - atoms[j].pos) < 2 * size and dot(atoms[i].pos - atoms[j].pos, atoms[i].v - atoms[j].v) < 0:
                atoms[i].v, atoms[j].v = vcollision(atoms[i], atoms[j])

    # Find collisions between the atoms and the walls, and handle their collisions
    for i in range(N):
        # Check collision with x walls
        if abs(atoms[i].pos.x) >= L_size:
            atoms[i].v.x = -atoms[i].v.x
            momentum += 2 * m * abs(atoms[i].v.x)

        # Check collision with y walls
        if abs(atoms[i].pos.y) >= L_size:
            atoms[i].v.y = -atoms[i].v.y
            momentum += 2 * m * abs(atoms[i].v.y)

        # Check collision with z walls
        if abs(atoms[i].pos.z) >= L_size:
            atoms[i].v.z = -atoms[i].v.z
            momentum += 2 * m * abs(atoms[i].v.z)

    if int(t / dt) % 1000 == 0:
        pressure = momentum / (6 * (2 * L) ** 2 * dt)/1000
        pressure_label.text = f'Pressure: {pressure:.2f} Pa'
        momentum = 0  # Reset momentum after each pressure calculation




20