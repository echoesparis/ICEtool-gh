# ICEtool-gh

Transposition of [ICEtool](https://github.com/Art-Ev/ICEtool) into a Grasshopper [Hops](https://github.com/mcneel/compute.rhino3d/tree/master/src/ghhops-server-py) python server for fast calculation of surface temperatures

We combined the initial 5-steps ICETool for QGIS workflow into a single function that you can deploy remotely with Flask. 

ICETool-gh provides a simplified component that allows to make informed design choices (e.g. vegetation, materials) that take into account Urban Heat Island (UHI) phenomena.



<p align="center">
<img src="images/ee_surface temperature_example.gif" title="example" />
</p>




Inputs : 

- Rhino .3dm file
- [Material database](https://github.com/Art-Ev/ICEtool/blob/main/Scripts/Example/Step_1/Material_database.csv) 

Hops function : 

* ee_surface-temperature.gh

Examples 

* ee_surface temperature_example.gh : use with ee_surface temperature_example.3dm
* ee_surface temperature_example_material database.gh

Outputs : 

- ground temperature.



<p align="center">
<img src="images/221209_02.jpg" title="gh-def" />
</p>



Ground temperature is an estimation based on :

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\pagecolor{white}Q_R=Q_H+Q_L+Q_C+{\delta}Q_S" title="ICEtool_computed" />
</p>

with:
- $Q_R$ : Heat flux related to radiation (from the sun, infrared radiation and the atmosphere)

- $Q_H$ : Heat flux related to convection (considered as very low and homogeneous)

- $Q_L$ : Sensitive and latent heat flux of water

- $Q_C$ : Heat flow related to conduction

- ${\delta}Q_S$ : Heat flow related to thermal storage (thermal capacity of materials)

  

### Install

With venv, pip, etc. like any other Python dev environment
```shell
    pip install requirements.txt
```
### Run

```shell
    python app.py
```



### Deploy

Follow [this tutorial](https://www.youtube.com/watch?v=SiCAIRc0pEI). The code is already set-up to work with Heroku.



### References

[What is Hops](https://developer.rhino3d.com/guides/compute/what-is-hops/)



