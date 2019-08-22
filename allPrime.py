from bokeh.models import widgets as wid
from bokeh.plotting import figure
from bokeh import plotting as plt
import bokeh.layouts as layout
from functools import partial
from bokeh.models import Legend
from numpy import linspace, asarray
from bokeh.models.glyphs import Patches, Text
from bokeh.models import ColumnDataSource

def findLens(data1, data2, numPrim):
    """
    Finds the lengths of all primer reads for the 2 coresponding experiments
    :param data1: the dictionary of all reads for the first experiment
    :param data2: the dictionary of all reads for the second experiment
    :param numPrim: the number of primer measured in both experiments.
    :return e1: a list with all lengths of all primers in the first experiment
    :return e2: a list with all lengths of all primers in the second experiment
    """
    e1 = []; e2 = []
    for p in range(numPrim):
        e1.append(len(data1['seq' + str(p)]) + 1) #add one just to make 0 counts show up nice on the log plot
        e2.append(len(data2['seq' + str(p)]) + 1)
    return(e1, e2)

def confInt(x, factor):
    """
    confInt will return a series of x and y points to define a polygon patch represting the confidence interval over
    the estimate for any given edge.
    """
    yBot = []; xBot = [] 
    yTop = []; xTop = [] #lists to put polygon intergers in
    for x_i in x:
        x_sd = (x_i)**(-1/2)
        y_sd = (factor*x_i)**(-1/2)
        yTop.append(factor*x_i + y_sd)
        xTop.append(x_i + x_sd)
        yBot.insert(0, factor*x_i - y_sd)
        xBot.insert(0, x_i - x_sd)
    
    return(yTop + yBot, xTop + xBot)

def allPrimePlot(e1, e2, e1Name, e2Name, names, diags, highlight):
    """
    Assembles the desired experiment data into a plot comparing relative abundance
    :param e1: the lengths of all primer reads in the 1st experiment
    :param e2: the lengths of all primer reads in the 2nd experiment
    :param e1Name: name of the first experiment
    :param e2Name: name of the second experiment
    :return p0: the Bokeh plot with desired features
    """
    p0 = figure(title = 'Comparison number of reads between experiments ' + e1Name + ' and ' + e2Name,
            plot_height = 700, plot_width = 800, tools = "hover, save, box_zoom, wheel_zoom, pan, reset",
            tooltips = [('Primer:', '@p'), ('Exp1 Reads: ', '@x'), ('Exp2 Reads: ', '@y')], 
            y_axis_type = "log", x_axis_type = "log", x_range = (.9*min(e1), 1.1*max(e1)),
            y_range = (.9*min(e2), 1.1*max(e2)))
    p0.xaxis.axis_label = e1Name
    p0.yaxis.axis_label = e2Name
    x = linspace(min(min(e1), min(e2)), max(max(e1), max(e2)), 80)
    legItems = []
    linColors = ['tomato', 'plum', 'salmon']
    for i, j in enumerate(diags):
        lin1Dict = dict(x = x, y = j*x)
        lin = p0.line(x = 'x', y = 'y', color = linColors[i%len(linColors)], line_dash = 'dashed', source = lin1Dict)
        legItems.append((str(j) + 'x', [lin]))
        xPatch, yPatch = confInt(x, j)
        c_source = ColumnDataSource(dict(x = xPatch, y = yPatch))
        cint = Patches(xs = 'x', ys = 'y')
        p0.add_glyph(c_source, cint)
    legend = Legend(items = legItems, location = (-80, 50))
    p0.add_layout(legend, 'right')
    p0.title.text_font_size = "15px"
    p0.title.align = 'center'
    col = []; alphas = []; textAlp = []; count = 0
    colors = ["cadetblue", "gold", "crimson", "cornflowerblue", "mediumorchid", "lightgreen", "dodgerblue", 
              "palevioletred", "yellowgreen", "hotpink"]
    for j in highlight:
        if j == 1:
            col.append(colors[count%len(colors)])
            alphas.append(0.7)
            textAlp.append(1)
            count = count + 1
        else:
            col.append("lightsteelblue")
            alphas.append(0.4)
            textAlp.append(0)
    cirDict = dict(p = names, x = e1, y = e2, col = col, alp = alphas)
    cir = p0.circle(x = 'x', y = 'y', size=20, color = "col", hover_color = 'aqua',
                    alpha = "alp", hover_alpha = 0.1, source = cirDict)
    lab_source = ColumnDataSource(dict(x = e1, y = e2, text = names, alpha = textAlp))
    labels = Text(x="x", y="y", text="text", text_alpha = "alpha")
    p0.add_glyph(lab_source, labels)
    return(p0)


def allPrimeTab(bDict, names, bNames):
    """
    When called from the main.py script, this assembles all of the componenets of the 
    """
    exp1 = 0; exp2 = 1 #default values at the beginning, we can change later.
    e1, e2 = findLens(bDict[exp1], bDict[exp2], len(names))
    diags = [1, 2, .5]
    highlight = [0 for i in range(len(names))]
    apPlot = allPrimePlot(e1, e2, bNames[exp1], bNames[exp2], names, diags, highlight)
    e1Drop = wid.Select(title = "Experiment 1:", value = bNames[0], options = bNames)
    e2Drop = wid.Select(title = "Experiment 2:", value = bNames[1], options = bNames)
    listBox = wid.MultiSelect(title = "Nodes to Highlight:", options = list(names), size = 20,
                              value = list(names))
    
    instructions = wid.Div(text = """<p style = "text-align:left;font-family:verdana"> 
        \t Select two experiments to compare. The plot \n
        will then display all of the reads across both \n
        primers between the 2 experiments selected.<p>""")
    
    #spacers to make the tab more visually appealing:
    h0 = wid.Div(height = 40); h1 = wid.Div(height = 50); h2 = wid.Div(height = 10); h3 = wid.Div(height = 20)
    h4 = wid.Div(height = 60); h5 = wid.Div(height = 10); h6 = wid.Div(height = 10); h7 = wid.Div(height = 20)
    w0 = wid.Div(width = 60); w1 = wid.Div(width = 30)
    
    replotButton = wid.Button(label = "Replot", button_type = "success")
    helpButton = wid.Button(label = "Help?")
    helpOptions = ["", "How do I read this plot?"]
    helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
    helpRead = wid.Div(width = 40)
    
    apTab = wid.Panel(child = layout.row(w0, layout.column(h0, instructions, h1, helpButton, h2, helpDrop, h3,
                                        helpRead), w1, apPlot, layout.column(h4, e1Drop, h5, e2Drop, h6, listBox,
                                                                             h7, replotButton)), 
                      title = "Compare Across Primers")
    passDict = {'tab': apTab, 'e1': e1Drop, 'e2': e2Drop, 'bDict': bDict, 'bNames': bNames, 'names': names, 'lb': listBox}
    replotButton.on_click(partial(allPrimeUpdate, apDict = passDict))
    helpDict = {'tab': apTab, 'selected': helpDrop, 'options': helpOptions}
    helpButton.on_click(partial(allPrimeHelp, hDict = helpDict))
    
    return(apTab)
    
def allPrimeUpdate(apDict):
    bNames = apDict['bNames']; bDict = apDict['bDict']; tab = apDict['tab']; names = apDict['names']
    exp1 = bNames.index(apDict['e1'].value); exp2 = bNames.index(apDict['e2'].value); on = apDict['lb'].value
    highlight = asarray([i in on for i in names])
    e1, e2 = findLens(bDict[exp1], bDict[exp2], len(names))
    diags = [1, 2, 0.5]
    updatePlot = allPrimePlot(e1, e2, bNames[exp1], bNames[exp2], names, diags, highlight.astype(int))
    tab.child.children[3] = updatePlot
     
def allPrimeHelp(hDict):
    
    tab = hDict['tab']; hReq = hDict['selected'].value; opts = hDict['options']
    
    #for each of the possible help requests, write a help message:
    if(hReq == opts[0]):
        tab.child.children[1].children[7] = wid.Div(text = """<p style = 
        text-align:left;font-family:verdana"> 
        \tPlease selecte a FAQ from the drop down to \n
        recieve help on!<p>""")
        
    if(hReq == opts[1]):
         tab.child.children[1].children[7] = wid.Div(text = """<p style = 
         "text-align:left;font-family:verdana"> 
         \t The plots are made on a log-log scale. This \n
         is just to make reading across multiple scales easier. \n
         The red line indicates where the points would have to \n
         lie to be at the same ratio, where are the two light \n
         blue lines represents 2x and 0.5x ratios respectively. \n
         \t In order to see which primers are outliers, move \n
         the mouse directly over them, and read the small box \n
         which comes up.<p>""")