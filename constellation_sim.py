import math
import os

import matplotlib.pyplot as plt
import numpy as np
# The path to the location of Basilisk
# Used to get the location of supporting data.
from Basilisk import __path__
from Basilisk.simulation import spacecraft
from Basilisk.utilities import SimulationBaseClass
from Basilisk.utilities import macros
from Basilisk.utilities import orbitalMotion
from Basilisk.utilities import simIncludeGravBody
from Basilisk.utilities import unitTestSupport  # general support file with common unit test functions
# attempt to import vizard
from Basilisk.utilities import vizSupport

bskPath = __path__[0]
fileName = os.path.basename(os.path.splitext(__file__)[0])

def check_target_range(chaser_num,target_num,orbit_para):
    if orbit_para[chaser_num-1]["i"] == orbit_para[target_num-1]["i"] and abs(orbit_para[chaser_num-1]["a"] - orbit_para[target_num-1]["a"])<5500000:
        chaser_a = orbit_para[chaser_num-1]["a"]
        target_a = orbit_para[target_num-1]["a"]
        print(f"chaser SMA : {chaser_a} and chaser satellite number is : {chaser_num}")
        print(f"targer SMA: {target_a} and target satellite number is : {target_num}")
        z =abs(orbit_para[chaser_num-1]["a"] - orbit_para[target_num-1]["a"])
        # print(f"here is the differnece in meet SMA{z}")
        print("you have enter the correct paramters the simulation will run sucessfully ")
        return 1
    else:
        print("chaser can't reach the target with the paramters priovided")
        return 0 


def run(show_plots, maneuverCase):
    # Create simulation variable names
    simTaskName = "simTask"
    simProcessName = "simProcess"

    #  Create a sim module as an empty container
    scSim = SimulationBaseClass.SimBaseClass()

    #
    #  create the simulation process
    #
    dynProcess = scSim.CreateNewProcess(simProcessName)

    # create the dynamics task and specify the integration update time
    simulationTimeStep = macros.sec2nano(10.)
    dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))
    scObjects = []
    orbit_para = {}
    position_list = []
    velocity_list = []
    rVt_list = []
    vVt_list = []

    while True:
        try:
            sat_num = int(input("Enter number of satellite in constellation [1-6] : ")) 
            if 1 <= sat_num <= 6:
                print(f"constellation is of : {sat_num} satellite ")
                break   # exit loop if valid
            else:
                print("Error: Number must be between 1 and 6.")
        except ValueError:
            print("Error: Please enter a valid integer.")

    gravFactory = simIncludeGravBody.gravBodyFactory()
    earth = gravFactory.createEarth()
    earth.isCentralBody = True  # ensure this is the central gravitational body   
    
# Create SPICE interface once
    timeInitString = "2024 September 21, 21:00:00.0 TDB"
    spiceObject = gravFactory.createSpiceInterface(time=timeInitString, epochInMsg=True)
    spiceObject.zeroBase = 'Earth'
    
    
# Add SPICE to sim
    scSim.AddModelToTask(simTaskName, spiceObject)


    for i in range(sat_num):
        sat = spacecraft.Spacecraft()  # create a new spacecraft instance
        sat.ModelTag = f"satellite{i+1}" # name/tag for the spacecraft
        scObjects.append(sat)           # store it in the list
        # add spacecraft object to the simulation process
        scSim.AddModelToTask(simTaskName, sat)
            

    # attach gravity model to spacecraft
        gravFactory.addBodiesTo(sat)
        oe = orbitalMotion.ClassicElements()
        print(f"\n--- Enter orbital parameters for Satellite {i+1} ---")
        a = float(input("Enter semi-major axis a (km): "))*1000
        e = float(input("Enter eccentricity e: "))
        i_deg = float(input("Enter inclination i (deg): "))
        Omega_deg = float(input("Enter RAAN Ω (deg): "))
        omega_deg = float(input("Enter argument of periapsis ω (deg): "))
        f_deg = float(input("Enter true anomaly f (deg): "))
        
        orbit_para[i]= {
            "a": a,
            "e": e,
            "i": i_deg * macros.D2R,
            "Omega": Omega_deg * macros.D2R,
            "omega": omega_deg * macros.D2R,
            "f": f_deg * macros.D2R
        }

    # Assign to Basilisk object
        oe.a = orbit_para[i]["a"]
        oe.e = orbit_para[i]["e"]
        oe.i = orbit_para[i]["i"]
        oe.Omega = orbit_para[i]["Omega"]
        oe.omega = orbit_para[i]["omega"]
        oe.f = orbit_para[i]["f"]

        rN, vN = orbitalMotion.elem2rv(earth.mu, oe)
        sat.hub.r_CN_NInit = rN
        sat.hub.v_CN_NInit = vN


            #
    #   Setup data logging before the simulation is initialized
    #
    # Compute simulation time from first satellite's orbit
    
    max_a = max(orbit_para[i]["a"] for i in range(sat_num))
    print("Maximum semi-major axis:", max_a)
    n = np.sqrt(earth.mu / max_a**3)
    P = 2. * np.pi / n
    print(f"Time period : {P}")
    simulationTime = macros.sec2nano(1.2*P)

    # Setup logging for all spacecraft
    for sat in scObjects:
        dataRec = sat.scStateOutMsg.recorder(simulationTimeStep)
        scSim.AddModelToTask(simTaskName, dataRec)

#choosing the satellites chaser parameters
    while True:
        try:
            chaser_num = int(input(f"Chose the chaser enter the satellite number (1 - {sat_num}): "))
            if 1 <=chaser_num <= sat_num:
                print(f"here is the number of meet : {sat_num}")
                print(f"your chaser satellite number is : {chaser_num}")
                break
            else:
                print(f"Error: Number must be between 1 and {chaser_num}.")
        except ValueError:
            print("Error: Please enter a valid integer.")
    #choosing the target satellite and checking targets parameters

    while True:
        try:
            target_num = int(input(f"Choose the target enter the satellite number other than satellite {chaser_num} : "))
            if 1<=target_num<= sat_num and target_num!=chaser_num:
                t = check_target_range(chaser_num,target_num,orbit_para)
                if t==1:
                    print(f"chaser satellite number is : {chaser_num} and target satellite number is : {target_num}")
                    break
                elif t==0:
                  check_target_range(chaser_num,target_num,orbit_para)
            # elif target_num == chaser_num:
            #     print("you have enter the chaser satellite")
            else:
                print("Error: Number must be between 1 and 6.amd it should not be same as chaser satellite")
        except ValueError:
            print("Error: Please enter valid target satellite parameters.")

    chaser = orbit_para[chaser_num-1]
    scObjects[chaser_num-1].ModelTag = "Chaser Satellite"
    scObjects[target_num-1].ModelTag = "Target Satellite"
    target = orbit_para[target_num-1]
    

    # Visualization (optional)
    if vizSupport.vizFound:
        vizSupport.enableUnityVisualization(
        scSim,
        simTaskName,
        scObjects,
        saveFile=fileName  # ✅ name of the file to save
    )
    # # ✅ Optional: direct line between satellites
    #     vizSupport.createLine(
    #        scSim,
    #        simTaskName,
    #        scObjects[chaser_num-1],
    #        scObjects[target_num-1],
    #        color=[0, 0, 1, 1],       # solid blue
    #        name="MWPT_Link"
    # )

    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    scSim.ExecuteSimulation()
    
    #for saving the data of position and velocity
    for sat in range(sat_num):
        posRef = scObjects[sat].dynManager.getStateObject(scObjects[sat].hub.nameOfHubPosition)
        velRef = scObjects[sat].dynManager.getStateObject(scObjects[sat].hub.nameOfHubVelocity)
        velocity_list.append(velRef)
        position_list.append(posRef)
        rVt_list.append(unitTestSupport.EigenVector3d2np(position_list[sat].getState()))
        print(f" here is the postion : {rVt_list[sat]}")
        vVt_list.append(unitTestSupport.EigenVector3d2np(velocity_list[sat].getState()))
        print(f" here is the velocity : {vVt_list[sat]}")
        print(f"Current Position : {rVt_list[sat]}")
        print(f"current Velocity : {vVt_list[sat]}")
        
    #checking wheather chaser is at perigee or not if yes then burn the thruster
    if orbit_para[chaser_num-1]["f"] == 0:
        v1=np.linalg.norm(vVt_list[chaser_num-1]) #chaser initial velocity
        r1=np.linalg.norm(rVt_list[chaser_num-1]) #chaser initial position
        r2=np.linalg.norm(rVt_list[target_num-1])#target initial position
        rhat_chaser = rVt_list[chaser_num-1]/r1#chaser direction 
        rhat_target = rVt_list[target_num-1]/r2#target direction
        hHat_chaser =np.cross(rVt_list[chaser_num-1],vVt_list[chaser_num-1])#chaser perpendicular direction
        hHat_target = np.cross(rVt_list[target_num-1],vVt_list[target_num-1])
        hHat_chaser= hHat_chaser/np.linalg.norm(hHat_chaser)
        hHat_target= hHat_target/np.linalg.norm(hHat_target)
        vHat_chaser = np.cross(hHat_chaser,rhat_chaser)#chaser velocity direction
        vHat_target = np.cross(hHat_target,rhat_target)
        at = (r1+r2)*0.5
        v1p=np.sqrt(earth.mu / r1 * r2 / at ) #tarsfer orbit velocity at chaser perigee
        vVt_list[chaser_num-1] = vVt_list[chaser_num-1] + vHat_chaser*(v1p - v1)
        n1 = np.sqrt(earth.mu/at/at/at)
        T2 = macros.sec2nano((np.pi) / n1)
        velocity_list[chaser_num-1].setState(vVt_list[chaser_num-1])
        orbit_para[chaser_num-1]["a"] = orbit_para[target_num-1]["a"] 
        orbit_para[chaser_num-1]["f"] = 180 


    scSim.ConfigureStopTime(simulationTime + T2)
    scSim.ExecuteSimulation()
# to save the current position after hohman transfer
    rVt_list[chaser_num-1] = (unitTestSupport.EigenVector3d2np(position_list[chaser_num-1].getState()))
    vVt_list[chaser_num-1] = (unitTestSupport.EigenVector3d2np(velocity_list[chaser_num-1].getState()))

    n2 = np.sqrt(earth.mu/r2/r2/r2)
    T3 = macros.sec2nano((np.pi) / n2)
    target_apogee = orbit_para[target_num-1]["a"]*(1 + orbit_para[target_num-1]["e"])
    error_in_apogee_burn = abs((np.linalg.norm(rVt_list[chaser_num-1])/1000) - (target_apogee/1000))
    if orbit_para[chaser_num-1]["f"] == 180 :
        v2 = np.linalg.norm(vVt_list[chaser_num-1]) #chaser initial velocity
        v2p = np.sqrt(earth.mu/r2) # velocity for the target orbit

        vHat_chaser = vVt_list[chaser_num-1]/v2
        vVt_list[chaser_num-1] = vVt_list[chaser_num-1]+ vHat_chaser*(v2p - v2)
        velocity_list[chaser_num-1].setState(vVt_list[chaser_num-1])
        
    print(f" here is the final position: {rVt_list[chaser_num-1]}")
    print(f" here is the final velocity: {vVt_list[chaser_num-1]}")
    r_vec = np.array(rVt_list[chaser_num-1])
    v_vec = np.array(vVt_list[chaser_num-1])
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    target_a = 1 / (2/r - v**2/earth.mu)
    print(f"Semi-major axis: {target_a/1000:.2f} km")
    chaser_a = orbit_para[chaser_num-1]["a"]
    error_in_SMA = 100*(chaser_a - target_a)/chaser_a
    print(f"here is the error in the SMA = {error_in_SMA}")


    scSim.ConfigureStopTime(simulationTime + T2 + 4*T3)
    scSim.ExecuteSimulation()




if __name__ == "__main__":
    run(
        True,  # show_plots
        0  # Maneuver Case (0 - Hohmann, 1 - Inclination)
    )
