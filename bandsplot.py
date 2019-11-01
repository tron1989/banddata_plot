import numpy as np
import matplotlib.pyplot as plt

hskllist=['A','B','C','D','E','F','G','H','I','J','K','L','M','N']

def plotbandone(filename='bandssoc.dat',ylims=(-1,1),efermi=0.0,ksg=7,kps=100,savefile='bands.eps',hsklname=[]):
    if hsklname==[]:    
        hsklname=hskllist[0:ksg+1]
    knum=ksg*kps
    bandsdata=np.loadtxt(filename)
    bandsnum=int(len(bandsdata)/knum)
    #calculate the band gap
    bandsp=[]  # positive
    bandsn=[]  # negative
    for value in bandsdata[:,1]:
        if value > 0:
            bandsp.append(value)
        else:
            bandsn.append(value)
    bandgap=np.array(bandsp).min()-np.array(bandsn).max()
    
    # hskl = high sysmestry kpoint location
    
    hskl=[bandsdata[0,0]]
    for value in range(ksg):
        hskl.append(bandsdata[(value+1)*kps-1,0])
    
    bandsper = np.zeros((knum,bandsnum))
    for i in range(bandsnum):
        bandsper[:,i]=bandsdata[i*knum:(i+1)*knum,1]-efermi
    
    plt.plot(bandsdata[0:knum,0],bandsper[:,[i for i in range(bandsnum)]],
    'b',linewidth='0.6')
    plt.ylim(ylims)
    plt.xlim((bandsdata[0:knum,0].min(),bandsdata[0:knum,0].max()))
    plt.axhline(y=0,linestyle='--',color='r',linewidth='0.6')
    for i in range(1,ksg+1):
        plt.axvline(x=hskl[i],linestyle='--',color='g',linewidth='0.6')
    plt.title('Bands Structure')
    plt.xlabel('gap: %.3f meV' %(bandgap*1000))
    plt.xticks(hskl,hsklname)
    plt.savefig(savefile,dpi=600,format='eps',transparent=True)
    plt.show()
        
def plotbandtwo(ksg=7,kps1=100,kps2=50,path1='bandnsoc.dat',path2='bandssoc.dat',ylims=(-2,2),savefile='bandswosoc.eps',hsklname=[]):
    
    # Working with multiple figure windows and subplots
    locatename=hsklname
    if locatename==[]:
        locatename=hskllist[0:ksg+1]
         
    pathnosoc=path1
    pathsoc=path2
#    ylims=(-2,2)
#    savefile='bandswosoc.eps'
#    #locatename=['Z','A','M',r'$\Gamma$','Z','R','X','V']
#    
#    #kpn=len(locatename)-1
    kpn=ksg
#    kps=100
    knumnosoc=kpn*kps1
    knumsoc=kpn*kps2
    
    bandsnosoc=np.loadtxt(pathnosoc)
    bandssoc=np.loadtxt(pathsoc)
    
    #for gap 
    socp=[]
    socn=[]
    
    for value in bandssoc[:,1]:
        if value>0:
            socp.append(value)
        else:
            socn.append(value)
    
    socgap=np.array(socp).min()-np.array(socn).max()
    
    nsocp=[]
    nsocn=[]
    
    for value in bandsnosoc[:,1]:
        if value>0:
            nsocp.append(value)
        else:
            nsocn.append(value)
    
    nsocgap=np.array(nsocp).min()-np.array(nsocn).max()
    
    
    
    #locatenosoc=[bandsnosoc[0,0],bandsnosoc[kps-1,0],bandsnosoc[2*kps-1,0],
    #             bandsnosoc[3*kps-1,0],bandsnosoc[4*kps-1,0]]
    
    locatenosoc=[bandsnosoc[0,0]]
    locatesoc=[bandssoc[0,0]]
    for i in range(len(locatename)-1):
        locatenosoc.append(bandsnosoc[((i+1)*kps1-1),0])
        locatesoc.append(bandssoc[((i+1)*kps2-1),0])
    #locatenosoc=[bandsnosoc[0,0],bandsnosoc[[((i+1)*kps-1) for i in range(len(locatename)-1)],0]]
    #locatesoc=[bandssoc[0,0],bandssoc[kps-1,0],bandssoc[2*kps-1,0],
    #             bandssoc[3*kps-1,0],bandssoc[4*kps-1,0]]
    
    
    knumnosocb=int(len(bandsnosoc)/knumnosoc)
    knumsocb=int(len(bandssoc)/knumsoc)
    
    bnosoc=np.zeros((knumnosoc,knumnosocb))
    for i in range(knumnosocb):
        bnosoc[:,i]=bandsnosoc[i*knumnosoc:(i+1)*knumnosoc,1]
        
    bsoc=np.zeros((knumsoc,knumsocb))
    for i in range(knumsocb):
        bsoc[:,i]=bandssoc[i*knumsoc:(i+1)*knumsoc,1]
        
    fig,(ax0,ax1)=plt.subplots(1,2,figsize=(15,5))
      
    #plt.figure(1)
    #plt.subplot(121)
    ax0.plot(bandsnosoc[0:knumnosoc,0], bnosoc[:,[i for i in range(knumnosocb)]],'b',linewidth='1.5')
    ax0.set_ylim(ylims)
    ax0.set_xlim((bandsnosoc[0:knumnosoc,0].min(),bandsnosoc[0:knumnosoc,0].max()))
    #plt.plot(bandsnosoc[0:knumnosoc,0], bnosoc[:,[i for i in range(knumnosocb)]],'b',linewidth='1.0')
    ax0.set_xticks(locatenosoc)
    ax0.set_xticklabels(locatename)
    ax0.axhline(y=0,linestyle='--',color='m',linewidth='0.6')
    for i in range(1,len(locatename)):
        ax0.axvline(x=locatenosoc[i],linestyle='--',color='g',linewidth='0.6')
    #    
    ax0.set_ylabel('Eenergy (eV)')
    ax0.set_xlabel('gap: %.3f meV' %(nsocgap*1000),fontdict={'fontweight' : 'bold'})
    ax0.set_title('Without SOC')
    #plt.subplot(122)
    #ax0.set_xlabel('(e) 4(8)-6-6-6-6')
    ax1.set_ylim(ylims)
    ax1.axhline(y=0,linestyle='--',color='r')
    ax1.set_xlim((bandssoc[0:knumsoc,0].min(),bandssoc[0:knumsoc,0].max()))
    ax1.plot(bandssoc[0:knumsoc,0], bsoc[:,[i for i in range(knumsocb)]],'r',linewidth='1.5')
    ax1.set_xticks(locatesoc)
    ax1.set_xticklabels(locatename)
    ax1.axhline(y=0,linestyle='--',color='m',linewidth='0.6')
    for i in range(1,len(locatename)):
        ax1.axvline(x=locatenosoc[i],linestyle='--',color='g',linewidth='0.6')
    ax1.set_xlabel('gap: %.3f meV' %(socgap*1000),fontname="Times New Roman")
    ax1.set_title('With SOC')
    #ax1.set_xlabel('(f) 4(8)-6-6-6-6')
    ax1.set_ylabel('Eenergy (eV)')
    #fig.set_label("Band structures")
    #
    fig.savefig(savefile,dpi=600,format='eps',transparent=True)
    plt.show()
    
def plotbandcomp(file1='bandnsoc.dat',file2='bandssoc.dat',
                 efermi1=0.0,efermi2=0.0,ksg=7,kwnk1=700,kwnk2=700,
                 ylims=(-1,1),savefile='bandscomp.eps',hsklname=[]):
    
#    file1='bandnsoc.dat'
#    file2='bandssoc.dat'
#    savefile='bandscomp.eps'
#    
#
#    efermi1=0.0
#    efermi2=0.0
#    ksg=7
#    kwnk1=700
#    kwnk2=700
#    ylims=(-1,1)
         
    if hsklname==[]:
        hsklname=hskllist[0:ksg+1]
        
    bands1=np.loadtxt(file1)
    bands2=np.loadtxt(file2)
    kwnk1a=int(len(bands1)/kwnk1)
    kwnk2a=int(len(bands2)/kwnk2)

    # hskl = high sysmestry kpoint location
    
    hskl=[bands1[0,0]]
    kps=int(kwnk1/ksg)
    for value in range(ksg):
        hskl.append(bands1[(value+1)*kps-1,0])
    

    # for the accidengt lines 
    bands1b=np.zeros((kwnk1,kwnk1a))
    for i in range(kwnk1a):
        bands1b[:,i]=bands1[i*kwnk1:(i+1)*kwnk1,1]-efermi1
        
    bands2b=np.zeros((kwnk2,kwnk2a))
    for i in range(kwnk2a):
        bands2b[:,i]=bands2[i*kwnk2:(i+1)*kwnk2,1]-efermi2
        
    
    #plt.plot(bands1[:,0],bands1[:,1])
    fig, ax1 = plt.subplots()
    #plt.ylim(-1,1)
    
    ax1.plot(bands1[0:kwnk1,0],bands1b[:,[i for i in range(kwnk1a)]],'b')
    ax1.set_xlim(bands1[:,0].min(),bands1[:,0].max())
    ax1.set_ylim(ylims)
    ax1.axhline(y=0,linestyle='--',color='r',linewidth='0.6')
    for i in range(1,ksg+1):
        ax1.axvline(x=hskl[i],linestyle='--',color='g',linewidth='0.6')
    ax1.set_xticks(hskl)
    ax1.set_xticklabels(hsklname)
    #ax1.set_xticks([bands1[0,0],bands1[99,0],bands1[199,0],bands1[299,0]])
    #ax1.set_xticklabels(hsklname)
#    ax1.axvline(x=bands1[99,0], color='g', linestyle='--' )
#    ax1.axvline(x=bands1[199,0], color='g', linestyle='--' )
    #for tl in ax1.get_xticklabels():
    #    tl.set_color('b')
    
    ax2=ax1.twiny()
    ax2.plot(bands2[0:kwnk2,0],bands2b[:,[i for i in range(kwnk2a)]],'r--')
    #ax2.plot(bands2[:,0],bands2[:,1:],'r--')
    ax2.set_ylim(ylims)
    ax2.set_xlim(bands2[:,0].min(),bands2[:,0].max())
    ax2.set_axis_off()
    #ax2.set_xticks([bands2[0,0],bands2[99,0],bands2[199,0],bands2[299,0]],
    #            [r'$M$',r'$\Gamma$',r'$K$',r'$M$'])
    
    plt.savefig(savefile,dpi=600,format='eps',transparent=True)
    
    plt.show()
    
def plotphonony(filepath='phonopyband.dat',ksg=3,kps=100,hsklname=[]):

    
    kwns=ksg
    kwnk=kps
    
    if hsklname==[]:
        hsklname=hskllist[0:ksg+1]
    
    file = open(filepath,'r')
    phonopy_bands=file.readlines()
    kwave=phonopy_bands[2].split()
    file.close()
    #for i in range(3):
    #    print(file.readline())
    #file.close()
    bands_data=np.loadtxt(filepath)
    
    bandlines=int(len(bands_data)/(kwns*kwnk))
    
    bandsb=np.zeros((kwns*kwnk, bandlines))
    for i in range(bandlines):
        bandsb[:,i]=bands_data[i*kwns*kwnk:(i+1)*kwns*kwnk,1]
    
    plt.plot(bands_data[0:kwns*kwnk,0],bandsb[:,[i for i in range(bandlines)]],'r',linewidth=3.0)
    #plt.axhline(y=0)
    plt.plot(bands_data[0:kwns*kwnk,0],np.zeros(kwns*kwnk),'k',linewidth=3.0)
    plt.xlim(bands_data[:,0].min(),bands_data[:,0].max())
    plt.ylim(0.0,bands_data[:,1].max())
    plt.xticks([float(kwave[i+1]) for i in range(len(kwave)-1)],
                hsklname)
          
    plt.savefig('phonony.eps',dpi=600,format='eps',transparent=True)
    plt.show()
