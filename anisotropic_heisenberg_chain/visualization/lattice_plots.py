import numpy as np
import plotly.graph_objects as go

NUM_ARROWS = 20
def create_1D_spin_tracer(ss: np.ndarray, color: str="blue", arrowhead: float=0.3,
                       opacity: float=0.5) -> list:
    maxlen = np.min([len(ss), NUM_ARROWS])    
    ss = ss[:maxlen]
    uu = np.arange(len(ss))
    vv = np.zeros(len(ss))
    ww = np.zeros(len(ss))
    xx = ss[:, 0]+uu
    yy = ss[:, 1]+vv
    zz = ss[:, 2]+ww
    x_lines = np.stack([uu,xx], axis=1)
    y_lines = np.stack([vv,yy], axis=1)
    z_lines = np.stack([ww,zz], axis=1)
    
    data = [
        
        go.Scatter3d(x=x_lines[i], y=y_lines[i], z=z_lines[i],
                    mode="lines", line=dict(color=color, width=6), opacity=opacity)
        
        for i in range(len(x_lines))
    ] +\
    [
        go.Cone(x=xx, y=yy, z=zz,
                u=ss[:,0], v=ss[:,1], w=ss[:,2],
                anchor="tip", showscale=False,
                colorscale=[[0, color], [1, color]],
                sizemode="scaled", sizeref=arrowhead, opacity=opacity),
        
        go.Scatter3d(x=uu, y=vv, z=ww, mode="markers",
                    marker=dict(size=10, color=color))
    ]
    return data

def make_1D_spin_plot(ss: np.ndarray, color: str="blue", arrowhead: float=0.3,
                        opacity: float=0.5) -> go.Figure:
        """Plot the spin configuration in 3D for each lattice site
    
        Args:
            ss (np.ndarray): The spin vectors at each lattice site (shape (N, 3))
            color (str, optional): Color for the vectors. Defaults to "blue".
            arrowhead (float, optional): Length of the arrow head. Defaults to 0.3.
            opacity (float, optional): Opacity of the objects. Defaults to 0.5.
    
        Returns:
            go.Figure: The plotly figure
        """
        data = create_1D_spin_tracer(ss, color=color,
                                     arrowhead=arrowhead, opacity=opacity)
        fig = go.Figure(data=data)
        fig.update_layout(scene=dict(xaxis=dict(range=[-1, len(ss)]),
                                    yaxis=dict(range=[-1, 1]),
                                    zaxis=dict(range=[-1, 1])),
                        width=700, margin=dict(r=20, l=10, b=10, t=10))
        fig.show()
        return None

def make_1D_spin_animation(ss_list: np.ndarray,  color: str="blue", arrowhead: float=0.3,
                      opacity: float=0.5) -> go.Figure:
    """Plot the spin configuration in 3D for each lattice site

    Args:
        ss (np.ndarray): The spin vectors at each lattice site (shape (N, 3))
        color (str, optional): Color for the vectors. Defaults to "blue".
        arrowhead (float, optional): Length of the arrow head. Defaults to 0.3.
        opacity (float, optional): Opacity of the objects. Defaults to 0.5.

    Returns:
        go.Figure: The plotly figure
    """

    n = len(ss_list)
    
    frames_data = [create_1D_spin_tracer(
        ss_list[idx], color=color,
        arrowhead=arrowhead, opacity=opacity) for idx in range(n)]
    
    step_list = np.arange(n)
    print("preparing animation...")
    fig = go.Figure(
        data=frames_data[0],
        layout=go.Layout(updatemenus=[dict(type="buttons",
                                           buttons=[dict(label="Play",
                                                         method="animate"
                                                         , args=[None])])], 
                         scene={'xaxis': {'range': [-1, NUM_ARROWS+1], 'rangemode': 'tozero', 'tickmode': "linear", 'tick0': -5, 'dtick': 1},
                                'yaxis': {'range': [-1, 1], 'rangemode': 'tozero', 'tickmode': "linear", 'tick0': -5, 'dtick': 1},
                                'zaxis': {'range': [-1, 1], 'rangemode': 'tozero', 'tickmode': "linear", 'tick0': -5, 'dtick': 1},
                            "aspectmode": "cube",
                             }),
                         
        frames=[go.Frame(data=frames_data[idx],
                         name="Step "+str(step_list[idx])) for idx in range(n)]
    )
    
    def frame_args(duration):
        return {
                "frame": {"duration": duration},
                "mode": "immediate",
                "fromcurrent": True,
                "transition": {"duration": duration, "easing": "linear"},
            }
    sliders = [
                {
                    "pad": {"b": 10, "t": 60},
                    "len": 0.9,
                    "x": 0.1,
                    "y": 0,
                    "steps": [
                        {
                            "args": [[f.name], frame_args(duration=0.1)],
                            "label": str(k),
                            "method": "animate",
                        }
                        for k, f in enumerate(fig.frames)
                    ],
                }
            ]
    fig.update_layout(sliders=sliders)
    fig.show()
    
    return None
