# Copyright 2013 Devsim LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from devsim.python_packages.simple_physics import *
from devsim.python_packages.ramp import *

import gmsh_mos2d_create

device = "mos2d"
silicon_regions=("gate", "bulk")
oxide_regions=("oxide",)
regions = ("gate", "bulk", "oxide")
interfaces = ("bulk_oxide", "gate_oxide")

set_parameter(device=device, name="initialTemperature", value=300.0)

for i in silicon_regions:
    SetSiliconParameters(device, i)

for i in oxide_regions:
    SetOxideParameters(device, i)


for i in regions:
    CreateSolution(device, i, "Potential")
    CreateSolution(device,i, "Temperature")
    CreateTinitial(device,i)
    set_node_values(device=device, region=i, name="Temperature", init_from="Init_temperature")
    edge_from_node_model(device=device, region=i, node_model="Temperature")
    edge_average_model(device=device, region=i, edge_model="edgeTemp", node_model="Temperature")
    edge_average_model(device=device, region=i, edge_model="edgeTemp", node_model="Temperature", derivative="Temperature")
##initial the temperature field for equilibrium solution

for i in silicon_regions:
    V_t_T = "boltzmannConstant*Temperature/(ElectronCharge)"
    node_model(device=device, region=i, name="V_t_T", equation=V_t_T)
    node_model(device=device, region=i, name="V_t_T:Temperature", equation="boltzmannConstant/ElectronCharge")
    edge_average_model(device=device, region=i, edge_model="edgeV_t_T", node_model="V_t_T")
    edge_average_model(device=device, region=i, edge_model="edgeV_t_T", node_model="V_t_T", derivative="Temperature")
##create model to calculate mobility and thermal voltage

for i in silicon_regions:
    CreateSiliconPotentialOnly(device, i)

for i in oxide_regions:
    CreateOxidePotentialOnly(device, i, "log_damp")

### Set up contacts
contacts = get_contact_list(device=device)
for i in contacts:
    tmp = get_region_list(device=device, contact=i)
    r = tmp[0]
    print("%s %s" % (r, i))
    CreateSiliconPotentialOnlyContact(device, r, i)
    # set_parameter(device=device, name=GetContactBiasName(i), value=0.0)
    set_parameter(device=device, name="drain_bias", value=0.0)
    set_parameter(device=device, name="source_bias", value=0.0)
    set_parameter(device=device, name="gate_bias", value=0.0)
    set_parameter(device=device, name="body_bias", value=0.0)

for i in interfaces:
    CreateSiliconOxideInterface(device, i)


solve(type="dc", absolute_error=1.0e-13, relative_error=1e-12, maximum_iterations=30)
solve(type="dc", absolute_error=1.0e-13, relative_error=1e-12, maximum_iterations=30)

write_devices(file="gmsh_mos2d_potentialonly", type="vtk")

#T field

for i in silicon_regions:
    CreateSolution(device, i, "Electrons")
    CreateSolution(device, i, "Holes")
    set_node_values(device=device, region=i, name="Electrons", init_from="IntrinsicElectrons")
    set_node_values(device=device, region=i, name="Holes",     init_from="IntrinsicHoles")
    # TODO
    set_parameter(device=device, region=i, name="mu_n", value=400.)
    set_parameter(device=device, region=i, name="mu_p", value=200.)
    set_parameter(device=device, region=i, name="T", value=300.)
    CreateMobilityModels(device, i, "mu_n", "mu_p")
    CreateSiliconDriftDiffusion(device, i)

#mu_n_T = "mu_n*pow((edgeTemp)/(2*T),-2.2)"
#mu_p_T = "mu_p*pow((edgeTemp)/(2*T),-2.2)"

for i in regions:
    ### TODO: create derivative at the same time as the temperature
    ### TODO: maybe edge_average_model doesn't need intermediate derivatives?
    if i in silicon_regions:
        CreateJouleHeatSilicon(device, i)
    elif i in oxide_regions:
        CreateJouleHeatOxide(device, i)
    else:
        raise RuntimeError("unknown region type for region " + i)

    CreateTemperaturefield(device,i)

for c in contacts:
    tmp = get_region_list(device=device, contact=c)
    r = tmp[0]
    CreateSiliconDriftDiffusionAtContact(device, r, c)
    # CreateTemperatureboundary(device, r, c)
    # set_parameter(device=device, name="drain_biasT", value=300)
    # set_parameter(device=device, name="source_biasT", value=300)
    # set_parameter(device=device, name="gate_biasT", value=300)
    # set_parameter(device=device, name="body_biasT", value=300)



CreateTemperatureboundary(device, "bulk", "body")  
set_parameter(device=device, name="body_biasT", value=300)  
# CreateTemperatureboundary(device, "bulk", "source")    
# set_parameter(device=device, name="source_biasT", value=300)
# CreateTemperatureboundary(device, "bulk", "drain")    
# set_parameter(device=device, name="drain_biasT", value=300)
# CreateTemperatureboundary(device, "gate", "gate")    
# set_parameter(device=device, name="gate_biasT", value=300)

##create heat sink
# CreateHeatfluxboundary(device,"gate","gate",0)
# CreateHeatfluxboundary(device,"bulk","source",0)
# CreateHeatfluxboundary(device,"bulk","drain",0)

# for c in ("source","drain"):
#     tmp = get_region_list(device=device, contact=c)
#     r = tmp[0]
#     CreateHeatfluxboundary(device,r,c,0)
##adiabatic boundary in contact

for i in interfaces:
    CreateTemperatureinterface(device, i)

solve(type="dc", absolute_error=1.0e30, relative_error=1e-5, maximum_iterations=1000)

for r in silicon_regions:
    node_model(device=device, region=r, name="logElectrons", equation="log(Electrons)/log(10)")

##write_devices -file gmsh_mos2d_dd.flps -type floops
##write_devices -file gmsh_mos2d_dd -type vtk
##write_devices -file gmsh_mos2d_dd.msh -type devsim_data

for r in silicon_regions:
    element_from_edge_model(edge_model="ElectricField",   device=device, region=r)
    element_from_edge_model(edge_model="ElectronCurrent", device=device, region=r)
    element_from_edge_model(edge_model="HoleCurrent",     device=device, region=r)


#
# printAllCurrents(device)
rampbias(device, "drain",  2, 0.01, 0.001, 100, 1e-10, 1e30, printAllCurrents)    ##current unit is A/cm
# rampbias(device, "gate", 1, 0.01, 0.001, 100, 1e-10, 1e30, printAllCurrents)
#
write_devices(file="gmsh_mos2d_dd.dat", type="tecplot")


