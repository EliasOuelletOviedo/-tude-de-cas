import numpy as np
import matplotlib.pyplot as plt

def graph(x,y, x_axis_name, y_axis_name):
    plt.plot(x,y, color="orange")
    plt.grid()
    plt.xlabel(x_axis_name, size=20)
    plt.ylabel(y_axis_name, size=20)
    plt.plot([0.75, 0.75], [min(y), max(y)], linestyle="dashed", color="gray")
    plt.plot([0.35, 0.35], [min(y), max(y)], linestyle="dashed", color="gray")

