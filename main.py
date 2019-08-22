from fileReading import *
from allPrime import *
from lenHist import *
from bokeh.models import widgets as wid
from bokeh.plotting import curdoc
from numpy import zeros, reshape, asmatrix, array
import pandas as pd
import bokeh.layouts as layout
from functools import partial
import os

def controlHelp(hDict):
    
    tab = hDict['tab']; hReq = hDict['selected'].value; opts = hDict['options']
    
    if(hReq == opts[0]):
        tab.child.children[7].children[4] = wid.Div(text = """<p style = 
        "text-align:center;font-family:verdana"> 
        Please selecte a FAQ from the drop down to recieve help on!<p>""")
        
    if(hReq == opts[1]):
         tab.child.children[7].children[4] = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
          This tab allows you to either run your own MCMC
          or load the results from a previous MCMC run for 
          analysis.  It defualts to loading an example file
          for you to play around with<p>""")
            
    if(hReq == opts[2]):
         tab.child.children[7].children[4] = wid.Div(text = """<p style = 
         "text-align:center;font-family:verdana"> 
          MCMC is a beautiful sampling algorithm.  When distributions
          are too hard to evaluate directly, either because they are 
          very high dimensional or the integral cannot be evaluated, 
          you have to approximate!  MCMC is a common way to approximate 
          these distribution (for example the one in this application)
          via sampling.  I encourage you to read more about it!<p>""")

def controlUpdate(pDict):
    gui = pDict['gui']; fileName = pDict['csv'].value; alpha = float(pDict['alpha'].value)
    burnIn = int(pDict['burnIn'].value); numSteps = int(pDict['numSteps'].value)
    thin = int(pDict['thin'].value); gamma = float(pDict['gamma'].value)
    
    dataDict, names = prepCSV(fileName)
    evolve, gs, lp, _, _ = minIMAP_MCMC(dataDict, burnIn = burnIn, alpha = alpha, justEdges = True, gamma = gamma, 
                 nSteps = numSteps, thinning = thin)
    edges = zeros(len(names)**2)
    for i in range(len(evolve[0,:])):
        cur = array(evolve[:,i]).flatten()
        edges[i] = cur[len(cur)-1]
    size = len(names)
    edges = reshape(asmatrix(edges), (size, size))
    
    gui.tabs[1] = fullTab(edges, names)    
    gui.tabs[2] = edgeTab(edges, names)
    gui.tabs[3] = learnTab(evolve, names)

cwd = os.getcwd()
primerFile = cwd + '\\primers.bowtie_txt_fixed'
alignDir = cwd + '\\Experiments\\190710_Rho_1.40_1.10_1.1\\aligned'
barcFile = cwd + '\\Experiments\\190710_Rho_1.40_1.10_1.1\\barcodeInfo.txt'

pDict = primerRead(primerFile)
bDict = processFiles(alignDir, pDict['loc'])
bNames = barcodeNames('', alignDir)

#sp0 = wid.Div(width = 60); sp1 = wid.Div(height = 40); 
#sp2 = wid.Div(width = 20); sp3 = wid.Div(height = 20) #spacers to make the tab more visually appealing
#sp4 = wid.Div(width = 40); sp5 = wid.Div(height = 40)
#sp6 = wid.Div(height = 100); sp7 = wid.Div(width = 40)
#sp8 = wid.Div(height = 40); sp9 = wid.Div(height = 60)
#sp10 = wid.Div(width = 40)
#
#instructions = wid.Div(text = """<p style = "text-align:left;font-family:verdana"> 
#    Welcome to the control panel! This panel allows you to navigate the entire 
#    process from raw allignment reads to data visualization.   If you have a 
#    specific question - see if it is in the FAQ dropdown menu and press 'Help?'
#    to ask that question! <p>""")
#sortStripText = wid.Div(text = """<p style = "text-align:left;font-family:verdana"> 
#    Sort and Strip the Poly A Tails from Raw Data: <p>""")
#alignText = wid.Div(text = """<p style = "text-align:left;font-family:verdana"> 
#    Align Sorted and Stripped Data: <p>""")
#plotText = wid.Div(text = """<p style = "text-align:left;font-family:verdana"> 
#    Plot Aligned Data: <p>""")
#
#alignButton = wid.Button(label = "Align", button_type = "success")
#sortStripButton = wid.Button(label = "Sort and Strip Raw Data", button_type = "success")
#plotButton = wid.Button(label = "Plot", button_type = "success")
#helpButton = wid.Button(label = "Help?", button_type = "success")
#
#primerFilInput = wid.TextInput(value = primerFile, title = "Location of Primer File:")
#alignDirInput = wid.TextInput(value = alignDir, title = "Directory of Aligned Files:")
#barcodeFile = wid.TextInput(value = barcFile, title = "Barcode Information File:")
#rawDataDir = wid.TextInput(value = '', title = "Raw Data Folder:")
#sortStripDir = wid.TextInput(value = '', title = "Sorted and Stripped Data Folder:")
#rawFileEnd = wid.TextInput(value = 'sequence.fastq', title = "Raw Data File Ending:")
#ssFileEnd = wid.TextInput(value = '.sorted.stripped', title = "Sorted and Stripped Data File Ending:")
#
#helpOptions = ["", "What should I do?", "What is MCMC?"]
#helpDrop = wid.Select(title = "FAQ:", value = helpOptions[0], options = helpOptions)
#helpRead = wid.Div(width = 40)

#tab0 = wid.Panel(child = layout.row(sp0, layout.column(sp1, sortStripText, sortStripDir, rawFileEnd,
#                                    ssFileEnd, sortStripButton), sp2, layout.column(sp3, alignText, 
#                                    sortStripDir, alignButton, plotText, alignDirInput, 
#                                    primerFilInput, plotButton), sp4,  layout.column(instructions, sp5,
#                                    helpButton, helpDrop, helpRead)),
#                title = "Control Panel")

#helpDict = {'tab': tab0, 'selected': helpDrop, 'options': helpOptions}
#helpButton.on_click(partial(controlHelp, hDict = helpDict))  

tab1 = allPrimeTab(bDict, pDict['names'], bNames)    
tab2 = lenHistTab(bDict, pDict['names'], bNames)
seqGUI = wid.Tabs(tabs = [tab1, tab2])

#passDict = {'gui': seqGUI, 'csv': csvInput, "alpha": alpha, 'burnIn': burnIn, 'numSteps': numSteps,
#           'thin': thinning, 'gamma': gamma}
#csvButton.on_click(partial(controlUpdate, pDict = passDict))

curdoc().add_root(seqGUI)
curdoc().title = "Target Seq Data Visualization"