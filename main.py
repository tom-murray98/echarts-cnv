import numpy as np
from nicegui import ui
from pathlib import Path
from cnv_from_bam import iterate_bam_file
import ruptures as rpt

# Assuming iterate_bam_file is defined elsewhere and works correctly

def show(event):
    ui.notify(f'{event.sender}: {event.value}')

# Path to BAM (currently hardcoded).
bam_path = Path("/home/mbxtm2/ds1305_Intraop0002_1.bam")

# Using cnv_from_bam to iterate over BAM and calculate CNV values per bin.
iterate_bam = iterate_bam_file(bam_path, _threads=4, mapq_filter=60)
bin_width = iterate_bam.bin_width

total = 0
penalty = 25
chr_list = []
scaled_changepoint_list = []


# Initialize an empty list to store all data series for plotting
series_data_list = []

for contig, cnv in iterate_bam.cnv.items():
    np_cnv = np.array(cnv)
    x_data = np.arange(len(cnv)) * bin_width + total
    y_data = cnv

    # Plotting the CNV data as scatter points or line
    series_data = {
        'name': f'{contig}',
        'type': 'scatter',  # Use 'scatter' for individual points, 'line' for continuous line
        'data': list(zip(x_data.tolist(), y_data)),
        'symbolSize': 3  # Adjust symbol size here
    }
    series_data_list.append(series_data)

    if len(np_cnv) > 3:
        # Change point detection
        algo = rpt.KernelCPD(kernel="linear", min_size=2).fit(np_cnv)
        result = algo.predict(pen=int(penalty))

        # Adding detected change points as red dashed lines
        for change_point in result:
            change_x = change_point * bin_width + total
            chr_list.append({'chromosome':contig,'breakpoint':change_x})
            series_data_list.append({
                'name': 'Change Point',
                'type': 'line',
                'data': [[change_x, min(y_data)], [change_x, max(y_data)]],
                'lineStyle': {'color': 'red', 'type': 'dashed', 'width': 2}
            })

    total += len(cnv) * bin_width
# Creating the plot in NiceGUI
ui.label('CNV Plot Showing Ruptures Breakpoints').classes('text-h5')
ui.echart({
    'tooltip': {'trigger': 'axis'},
    'xAxis': {'type': 'value'},
    'yAxis': {'type': 'value', 'max': 5, 'min': -1},
    'dataZoom': [  # Enable data zoom component
        {'type': 'inside'},  # Enable zooming using mouse wheel
        {'type': 'slider'}   # Enable zooming by selecting a region
    ],

    'series': series_data_list
})


columns = [
    {'name': 'chromosome', 'label': 'Chromosome', 'field': 'chromosome', 'required': True, 'align': 'left'},
    {'name': 'breakpoint', 'label': 'Breakpoint', 'field': 'breakpoint', 'sortable': True},
]
rows = chr_list

ui.table(columns=columns, rows=rows, row_key='name', pagination=5)

ui.run(port=8895, host='0.0.0.0')

