from dash import Dash, html, dcc, Input, Output  # pip install dash
import dash_ag_grid as dag
import dash_bootstrap_components as dbc   # pip install dash-bootstrap-components

import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO

import planetary_data
from plot_orbits_interface import plot_orbits
from tools import get_satellite_names, create_dataframe_and_write_to_csv

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Example usage
folder_path = '../data'  # Replace with the actual path to your data folder

satellite_names = get_satellite_names('../data')
satnames = create_dataframe_and_write_to_csv(satellite_names).sort_values(by='Satellite Names')
print(satnames)

app.layout = dbc.Container([
    html.H1("Interactive Matplotlib with Dash", className='mb-2', style={'textAlign':'center'}),
    dcc.Dropdown(
                satnames['Satellite Names'],
                'CHINASAT 18',
                id='active_sat',
                style={"width": "50%"},
            ),
    dbc.Row([
        dbc.Col([
            html.Img(id='bar-graph-matplotlib')
        ], width=12)
    ]),

])

# Create interactivity between dropdown component and graph
@app.callback(
    Output(component_id='bar-graph-matplotlib', component_property='src'),
    Input('active_sat', 'value'),
)
def plot_data(sat):
    # sats=['CHINASAT 18','YAOGAN-40 01C','BEIDOU 3 G4']
    sats = [sat]
    # perturbations=[[planetary_data.mars]]
    perturbations = ['j2', [planetary_data.moon]]
    bodies = []
    # bodies=['earth', 'mars']
    t0 = 'Jul 31, 2023, 00:00 UTC'
    tf = 'Aug 07, 2023, 00:00 UTC'
    cb = planetary_data.earth

    sc_data = {
        'mass': 5000,  # kg
        'reflectance': 0.2,
        'area': 12  # m²
    }

    fig = plot_orbits(sats, bodies, t0, tf, perturbations, central_body=cb, sc_data=sc_data, steps=4000, animate=False,
                show=True)


    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,{fig_data}'


    return fig_bar_matplotlib


if __name__ == '__main__':
    app.run_server(debug=True)