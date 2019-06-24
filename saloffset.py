import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from functools import partial
import nstools

def submit(profiles,cruisename,fig,ax,refindex,indexs,s,text):
    val = eval(text)
    ax.clear()
    ax.plot(profiles[refindex].isals,profiles[refindex].itemps,color="red",linewidth=3.0)
    if indexs == []:
        for p in range(len(profiles)):
            if profiles[p].cruise == cruisename:
                ax.plot(profiles[p].isals+val,profiles[p].itemps)
                indexs.append(p)
    else:
        for i in indexs:
            if profiles[i].cruise == cruisename:
                ax.plot(profiles[i].isals+val,profiles[i].itemps)

    ax.set_xlabel("salinity") 
    ax.set_ylabel("temperature") 
    plt.draw()
    s.answer = val
    return indexs

def selectorGraph(profiles,cruisename,refindex):
    fig, ax1 = plt.subplots(1,1)
    plt.subplots_adjust(bottom=0.2)
    indexs = submit(profiles,cruisename,fig,ax1,refindex,[],selectorGraph,"0.0")
    axbox = plt.axes([0.1, 0.05, 0.2, 0.04])
    text_box = TextBox(axbox, 'Evaluate', initial="0.0")
    text_box.on_submit(partial(submit,profiles,cruisename,fig,ax1,refindex,indexs,selectorGraph))
    plt.show()
    return selectorGraph.answer
    


profiles,deepestindex = nstools.extractProfilesMonths('data/2008profiles.json',range(13))

print(selectorGraph(profiles,profiles[0].cruise,0))

#print("DONE WITH EXTRACTING PROFILES")
#surfaces = search(profiles,deepestindex)
#print("DONE FINDING SURFACES")
#with open('data/surfaces.json', 'w') as outfile:
    #json.dump(surfaces, outfile)
##json_file = open("data/surfaces.json") 
##surfaces = json.load(json_file)
#print("NOW GRAPHING")
#graphSurfaces(surfaces)
