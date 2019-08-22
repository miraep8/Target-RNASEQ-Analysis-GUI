from bokeh.models import widgets as wid
from bokeh.plotting import figure
from bokeh import plotting as plt
import bokeh.layouts as layout
from functools import partial
from bokeh.models import Legend
from numpy import histogram, zeros


def makeHist(dat):
    """
    given the data, make a histogram to plot
    :param dat: the data to transform into a histogram
    :return hist: the heights of all histogram rects
    :return edge: the edges of all histogram rects
    """
    x = []
    for i in dat:
        if len(i) != 50:
            x.append(len(i))
    hist, edge = histogram(x)
    return(hist, edge)

def lenHistPlot(primName, expName, hist, edges):
    """
    Assembles the desired experiment data into a plot comparing relative abundance
    :param primName: the name of the primer currently being plotted
    :param expName: the name of the current experiment
    :param hist: the heights of the histograms
    :param edges: the edges of the histogram
    :return p: the Bokeh plot with desired features
    """
    p = figure(title = 'Lengths of reads for primer ' + primName + ' for Barcode ' + expName,
            plot_height = 600, plot_width = 600, tools = "hover, save, box_zoom, wheel_zoom, pan, reset",
            tooltips = [('Counts:', '@t'), ('Length:', '@l' ' - ' '@r')])
    hDict = dict(t = hist, b = zeros(hist.shape), l = edges[:-1], r = edges[1:])
    h = p.quad(top = 't', bottom = 'b', left = 'l', right = 'r', color = '#41b6c4', line_color = '#225ea8',
                 alpha = 0.7, hover_alpha = 0.8, source = hDict)
    p.xaxis.axis_label = 'Length of Read'
    p.yaxis.axis_label = 'Number of Reads'
    p.title.text_font_size = "15px"
    p.title.align = 'center'
    return(p)

def lenHistTab(bDict, names, bNames):
    """
    When called from the main.py script, this assembles all of the componenets of the 
    """
    exp1 = 0; exp2 = 1; dat1 = bDict[exp1]; dat2 = bDict[exp2] #default values at the beginning, we can change later.
    h1, e1 = makeHist(dat1['seq0']); h2, e2 = makeHist(dat2['seq0'])
    hist1 = lenHistPlot(names[0], bNames[0], h1, e1)
    hist2 = lenHistPlot(names[0], bNames[1], h2, e2)
    e1Drop = wid.Select(title = "Experiment 1:", value = bNames[0], options = bNames)
    e2Drop = wid.Select(title = "Experiment 2:", value = bNames[1], options = bNames)
    pDrop = wid.Select(title = 'Primer:', value = names[0], options = names)
    
    instructions = wid.Div(text = """<p style = "text-align:left;font-family:verdana"> 
        \t Select two experiments to compare from below \n 
        the plots. Then to the left of the plots select \n
        the primer you wish to make a comparison on.<p>""")
    
    #spacers to make the tab more visually appealing:
    h0 = wid.Div(height = 40); h1 = wid.Div(height = 50); h2 = wid.Div(height = 40); h3 = wid.Div(height = 10)
    h4 = wid.Div(height = 5); h5 = wid.Div(height = 10); 
    w0 = wid.Div(width = 60)
    
    replotButton = wid.Button(label = "Replot", button_type = "success")
    helpButton = wid.Button(label = "Help?")
    helpOptions = ["", "How do I read this plot?"]
    helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
    helpRead = wid.Div(width = 40)
    
    lhTab = wid.Panel(child = layout.row(w0, layout.column(layout.row(hist1, hist2), 
                                        layout.row(e1Drop, e2Drop, pDrop, replotButton)), 
                                        layout.column(h2, instructions, h3, helpDrop, h4, helpButton, h5,
                                        helpRead)), 
                      title = "Length Distribution Comp")
    passDict = {'tab': lhTab, 'e1': e1Drop, 'e2': e2Drop, 'p': pDrop, 'bDict': bDict, 
                'bNames': bNames, 'names': names}
    replotButton.on_click(partial(lenHistUpdate, lhDict = passDict))
    helpDict = {'tab': lhTab, 'selected': helpDrop, 'options': helpOptions}
    helpButton.on_click(partial(lenHistHelp, hDict = helpDict))
    
    return(lhTab)
    
def lenHistUpdate(lhDict):
    bNames = lhDict['bNames']; bDict = lhDict['bDict']; tab = lhDict['tab']; names = lhDict['names']
    exp1 = bNames.index(lhDict['e1'].value); exp2 = bNames.index(lhDict['e2'].value);
    p = names.index(lhDict['p'].value)
    
    dat1 = bDict[exp1]; dat2 = bDict[exp2] #default values at the beginning, we can change later.
    h1, e1 = makeHist(dat1['seq' + str(p)]); h2, e2 = makeHist(dat2['seq' + str(p)])
    hist1 = lenHistPlot(names[p], bNames[exp1], h1, e1)
    hist2 = lenHistPlot(names[p], bNames[exp2], h2, e2)
    
    tab.child.children[1].children[0].children[0] = hist1
    tab.child.children[1].children[0].children[1] = hist2
     
def lenHistHelp(hDict):
    
    tab = hDict['tab']; hReq = hDict['selected'].value; opts = hDict['options']
    
    #for each of the possible help requests, write a help message:
    if(hReq == opts[0]):
        tab.child.children[2].children[7] = wid.Div(text = """<p style = 
        text-align:left;font-family:verdana"> 
        \tPlease selecte a FAQ from the drop down to \n
        recieve help on!<p>""")
        
    if(hReq == opts[1]):
         tab.child.children[2].children[7] = wid.Div(text = """<p style = 
         "text-align:left;font-family:verdana"> 
         \t To Do!<p>""")