# Presentation

A Grasshopper [Hops](https://github.com/mcneel/compute.rhino3d/tree/master/src/ghhops-server-py) python server that easily computes ground temperatures in an urban environment.

This allows you to make informed design choices (e.g. vegetation, materials) that reduce urban heat island phenomena.

<p align="center">
<img src="images/ee_surface temperature_example.gif" title="example" />
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

[ICEtool](https://github.com/Art-Ev/ICEtool)

[What is Hops](https://developer.rhino3d.com/guides/compute/what-is-hops/)



