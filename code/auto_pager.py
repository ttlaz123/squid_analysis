import os 

def script_pad():
    script_str = "function pad(num, size) {\n"
    script_str += "    num = num.toString();\n"
    script_str += "    while (num.length < size) num = \"0\" + num;\n"
    script_str += "    return num;\n"
    script_str += "}\n\n"
    return script_str


def script_rows(numrows):
    script_str = "const rows = [];\n"
    script_str += f"numrows = {numrows}\n"
    script_str += "for (let i = 0; i < numrows; i++) {\n"
    script_str += "    s = pad(i, 1);\n"
    script_str += "    rows[i] = s.concat('|', s);\n"
    script_str += "}\n\n"
    script_str += "rows[numrows] = 'Summary';\n\n"
    return script_str


def script_cols(startcols, numcols):
    script_str = "const cols = [];\n"
    script_str += f"startcols = {startcols}\n"
    script_str += f"numcols = {numcols}\n"
    script_str += "for (let i = startcols; i < numcols; i++) {\n"
    script_str += "    s = pad(i, 1);\n"
    script_str += "    cols[i - startcols] = s.concat('|', s);\n"
    script_str += "}\n\n"
    return script_str


def script_link(ctime, figs_dir):
    script_str = "pager.link(\"#ic\",\n"
    script_str += "    {\n"
    script_str += "        'Col|col': cols,\n"
    script_str += "        'Row|row': rows,\n"
    script_str += f"        'Run|run': [\n"
    script_str += f"            '{ctime}|{ctime}',\n"
    script_str += "        ],\n"
    script_str += "        'Units|units': [\n"
    script_str += "            'Microamps (warning, calibration very likely incorrect)|uA',\n"
    script_str += "            'DAC|dac'\n"
    script_str += "        ]\n"
    script_str += "    },\n"
    script_str += "    function (params) {\n"
    script_str += "        ctime = params.run;\n"
    script_str += "        if (params.row == 'Summary') {\n"
    script_str += "            row = '_summary';\n"
    script_str += "            row_folder = 'col_summary';\n"
    script_str += "        } else {\n"
    script_str += "            row = '_row' + params.row;\n"
    script_str += "            row_folder = 'all_rows';\n"
    script_str += "        }\n"
    script_str += "        units = params.units;\n"
    script_str += "        if (units == 'uA') {\n"
    script_str += "            units = 'ua';\n"
    script_str += "            units_folder = 'units_ua';\n"
    script_str += "        } else if(units == 'dac'){\n"
    script_str += "            units = 'DAC';\n"
    script_str += "            units_folder = 'units_dac';\n"
    script_str += "        }\n"
    script_str += f"        name = \"{figs_dir}\" + ctime + \"/\" + units_folder + \"/\" + row_folder + \"/\" + ctime + '_icminmax_units' + units + '_col' +\n"
    script_str += "            params.col  + row+ '.png';\n"
    script_str += "        console.log(name);\n"
    script_str += "        return name;\n"
    script_str += "});\n"
    return script_str

def gridplots(ctime, figs_dir):
    script_str = "<figure>\n"
    script_str += "    <img alt=\"Squid Tuning Paramter Grid Plots\" id=\"grid\" src=\"#\" onerror=\"javascript:this.src='dne.png'\" />\n"
    script_str += " <figcaption>\n"
    script_str += "     <p>\n"
    script_str += "Grid plots\n"
    script_str += "     </p>\n"
    script_str += " </figcaption>\n"
    script_str += "<script type=\"text/javascript\">\n\n"
    script_str += f"pager.link(\"#grid\",\n"
    script_str += "    {\n"
    script_str += "        'Grid Type|grid':[\n"
    script_str += "            'Modulation at Chosen Bias|chosenmod',\n"
    script_str += "            'Ic,col|ic_col',\n"
    script_str += "            'Ic,max|ic_max',\n"
    script_str += "            'Optimal Bias|optbias',\n"
    script_str += "            'Modulation at optimal bias|optmod',\n"
    script_str += "            'Crosstalk Bias Limit|crosstalk',\n"
    script_str += "            'Ic,col-Ic,max|ic_maxcol_diff',\n"
    script_str += "            'Crosstalk bias limit - optimal bias|optbias_crosstalk_diff',\n"
    script_str += "        ],\n"
    script_str += "        'Run|run': [\n"
    script_str += f"            '{ctime}|{ctime}',\n"
    script_str += "        ],\n"
    script_str += "        'Units|units': [\n"
    script_str += "            'Microamps (warning, calibration very likely incorrect)|uA',\n"
    script_str += "            'DAC|dac'\n"
    script_str += "        ]\n"
    script_str += "    },\n"
    script_str += "    function (params) {\n"
    script_str += "        ctime = params.run;\n\n"
    script_str += "        row_folder = 'gridplots';\n\n"
    script_str += "        units = params.units;\n"
    script_str += "        if (units == 'uA') {\n"
    script_str += "            units = 'ua';\n"
    script_str += "            units_folder = 'units_ua';\n"
    script_str += "        } else if(units == 'dac'){\n"
    script_str += "            units = 'unitsDAC';\n"
    script_str += "            units_folder = 'units_dac';\n"
    script_str += "        }\n"
    script_str += f"        name = \"{figs_dir}/\" + ctime + \"/\" + units_folder + \"/\" + row_folder + \"/\" + ctime + \n"
    script_str += "                '_' + params.grid + '_' + units +'.png';\n"
    script_str += "        console.log(name);\n"
    script_str += "        return name;\n"
    script_str += "});\n"
    script_str += f"pager.setparams(\n"
    script_str += "    {\n"
    script_str += "        'grid': 'ic_maxcol_diff',\n"
    script_str += f"        'run': '{ctime}',\n"
    script_str += "        'units': 'dac'\n"
    script_str += "    }\n"
    script_str += ");\n"
    script_str += "    </script>\n"
    script_str += "</figure>"
    
    return script_str

def script_setparams(ctime):
    script_str = "pager.setparams(\n"
    script_str += "    {\n"
    script_str += "        'col': '0',\n"
    script_str += "        'row': 'Summary',\n"
    script_str += f"        'run': '{ctime}',\n"
    script_str += "        'units': 'dac'\n"
    script_str += "    }\n"
    script_str += ");\n\n"
    return script_str


def script(ctime, rows, cols, figs_dir):
    
    if(isinstance(rows, int)):
        numrows = rows
    else:
        numrows = len(rows)
    if(isinstance(cols, int)):
        numcols = cols 
        startcols = 0
    else:
        numcols = len(cols)
        startcols = min(cols)
    script_str = "<script type=\"text/javascript\">\n"
    script_str += script_pad()
    script_str += script_rows(numrows)
    script_str += script_cols(startcols, numcols)
    script_str += script_link(ctime, figs_dir)
    script_str += script_setparams(ctime)
    script_str += "</script>\n"
    return script_str


def caption(caption=None):
    if(caption is None):
        caption =("Green/aqua lines: Row selects off; when the green and aqua lines separate, that"+
                    "indicates that the columns is turning on even though the switches are off."+
                    "Red/blue lines: row selects on. Dark purple vertical line"+
                    "represents the amplitude of the maximum modulation when row select is on."+
                
                "The dotted thick pink line is the"+
                    "crosstalk threshold and marks the bias current and device current where"+
                    "a column turns on even when all the row selects are off. This is a signature of"+
                    "potential crosstalk. The small pink line represents the manually chosen bias"+
                    "current to maximize the modulation while keeping the SQ1 Max below the"+
                    "crosstalk threshold (small pink line only exists on summary plots).")
    cap_str = " <figcaption>\n"
    cap_str += "     <p>\n" + str(caption) + "\n      </p>\n" 
    cap_str += " </figcaption>\n"

    return cap_str

def auto_pager(rows, cols, ctime, convert_units=False, 
               figs_dir=None, html_dir='webpage'):
    '''
    Writes a pager for viewing the plots
    Requires the pager.css and pager.js files in the repo
    TODO In progress
    '''
    if(figs_dir is None):
        figs_dir = '../output_data/'
    
    figure_str = "<!DOCTYPE html>\n"
    figure_str += "<body>\n"
    figure_str +='<script type="text/javascript" src="scripts/pager.js"></script>\n'
    figure_str +='<link rel="stylesheet" type="text/css" href="scripts/pager.css">\n'
    figure_str += "<figure>\n"
    figure_str += (" <img alt=\"Squid Tuning Paramters\" "+
                   "id=\"ic\" src=\"#\" onerror=\"javascript:this.src='dne.png'\" />\n")
    figure_str += caption()

    figure_str += script(ctime, rows, cols, figs_dir)
    figure_str += "</figure>\n"
    figure_str += gridplots(ctime, figs_dir)
    figure_str += "</body>\n"

    os.makedirs(html_dir, exist_ok=True)

    html_filename = 'index.html'
    html_file_path = os.path.join(html_dir, html_filename)
    with open(html_file_path, 'w') as html_file:
        html_file.write(figure_str)


    return 

def main():
    ## test run

    auto_pager(rows=41, cols=16, ctime=1684269423)

if __name__ == '__main__':
    main()