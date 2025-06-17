import streamlit as st
import numpy as np
import math
import matplotlib.pyplot as plt
from orbittools import *

st.set_page_config(layout="wide")

col1, col2 = st.columns([1, 2])
with col1:
    st.subheader("Inputs:", divider="blue")

    c1, c2, c3 = st.columns(3)
    with c1:
        in1 = st.number_input("Position 1 (AU):", value=1.0, format="%0.3f", key=1)
    with c2:
        in2 = st.number_input("Y:", label_visibility="hidden", value=0.0, format="%0.3f", key=2)
    with c3:
        in3 = st.number_input("Z:", label_visibility="hidden", value=0.0, format="%0.3f", key=3, disabled=True)

    c4, c5, c6 = st.columns(3)
    with c4:
        in4 = st.number_input("Position 2 (AU):", value=0.0, format="%0.3f", key=4)
    with c5:
        in5 = st.number_input("Y:", label_visibility="hidden", value=1.0, format="%0.3f", key=5)
    with c6:
        in6 = st.number_input("Z:", label_visibility="hidden", value=0.0, format="%0.3f", key=6, disabled=True)

    mu_def = math.pi**2 * 4
    mu = mu_def

    r1 = np.array([in1, in2, in3])
    r2 = np.array([in4, in5, in6])
    if np.linalg.norm(np.cross(r1, r2)) == 0:
        hohmann = True
        a = (np.linalg.norm(r1)+np.linalg.norm(r2))/2
        tof_def = 0.0
    else:
        hohmann = False
        tof_def = 0.25
    
    tof = st.number_input("Time of Flight (yrs):", 0.0, key=7, value=tof_def, disabled=hohmann)
    mu = st.number_input("Gravitational Parameter:", 0.01, key=8, value=mu_def)

    plottype = st.checkbox("Only show transfer trajectory")

v1, v2 = lambert(r1, r2, mu, tof)
hv = np.cross(r1, v1)
E = np.linalg.norm(v1)**2/2-mu/np.linalg.norm(r1)
a = -mu/(2*E)
ev = np.cross(v1, hv)/mu - r1/np.linalg.norm(r1)
if np.linalg.norm(np.cross(ev, r1)) == 0:
    p = np.linalg.norm(hv)**2/mu
    e = math.sqrt(1-p/a)
    if np.linalg.norm(r1) < a:
        periapsis = r1
        w = 0
    else:
        periapsis = -r1 * (a*(1-e)/np.linalg.norm(r1))
        w = math.pi
else:
    e = np.linalg.norm(ev)
    p = a*(1-e**2)
    periapsis = ev/e * (p/(1+e))
    w = math.acos(periapsis[0]/np.linalg.norm(periapsis)) * (-periapsis[1]/abs(periapsis[1]))
if np.linalg.norm(ev) < 1e-12:
        w = 0
        e = 0

with col2:
    fig = plt.figure()
    if plottype:
        plotSolution(r1, r2, v1, mu, 10000, (1, 0, 1))
    else:
        plotOrbit(r1, v1, mu, 10000, (1, 0, 1))
    plt.plot(np.array([0, r1[0]]), np.array([0, r1[1]]))
    plt.plot(np.array([0, r2[0]]), np.array([0, r2[1]]))
    plt.plot(0, 0, 'xk')
    plt.xlabel("X (AU)")
    plt.ylabel("Y (AU)")
    st.pyplot(fig)

with col1:
    st.subheader("Orbital Elements:", divider="blue")
    st.write("Semimajor axis: ", a)
    st.write("Eccentricity: ", e)
    st.write("Angle from Periapsis: ", w)
