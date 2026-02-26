# Satellite Constellation Rendezvous Simulation
This project simulates a satellite constellation and performs feasibility analysis for orbital rendezvous between two selected spacecraft. The simulation is built using Basilisk for astrodynamics modeling and Vizard for 3D visualization.

The user can dynamically define the constellation size, orbital parameters, and select a chaser–target pair. If the rendezvous is feasible based on orbital conditions, the simulation runs and visualizes the scenario in Vizard.

# **Requirements**

Before running the project, ensure the following are installed:

1. Basilisk
2. Vizard (Basilisk Visualization Tool)
3. Python (compatible with your Basilisk installation)
You must properly source the Basilisk environment before execution.

# **Activate Basilisk Environment**

**Windows (PowerShell):**    

.\venv\Scripts\activate

**Linux / macOS:**

source venv/bin/activate

# **How to Run the Simulation**

py constellation.py

# **Simulation Workflow**

1. The program first prompts you to enter:

2. Number of satellites in the constellation.

3. For each satellite, you will provide:

4. Orbital parameters (e.g., semi-major axis, eccentricity, inclination, etc.).

5. After defining the constellation:

6. Select the chaser satellite (by index number).

7. Select the target satellite (by index number).

8. The system checks Whether orbital transfer/rendezvous is feasible.

9. If feasible the satellite will perform the orbit change and if not then need to re-enter the orbital parameter.

10. You can observe the 3D simulation of the rendezvous using loading .bin file in Vizard.
