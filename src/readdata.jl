    using PyPlot

    x = admmst.x
    mydata = skyst.mydata
    sky = skyst.sky
    nu = psfst.nu
    nu0 = psfst.nu0

    #############################################################################
    #############################################################################
    ######################################  sky / dirty / rec_obj / error
    figure(1)
    clf()
    title("multiplot")
    e = 1
    for z in [1 5 10 15]
        subplot(4,4,e)
        axis("off")
        a = trunc(nu[z]/1e9,3)
        imshow(sky[:,:,z],vmin=0,vmax=1.3)
        title(string(L"$\nu \, = \,$","$a ", L"\, GHz"),fontsize=14)
        e += 1
    end
    e = 5
    for z in [1 5 10 15]
        subplot(4,4,e)
        axis("off")
        imshow(mydata[:,:,z]./maximum(mydata).*1.2,vmin=0,vmax=1.3)
        e += 1
    end
    e = 9
    for z in [1 5 10 15]
        subplot(4,4,e)
        axis("off")
        imshow(x[:,:,z],vmin=0,vmax=1.3)
        e += 1
    end

    subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = axes([0.85, 0.3, 0.025, 0.6])
    colorbar(cax=cax)


    e = 13
    for z in [1 5 10 15]
        subplot(4,4,e)
        axis("off")
        imshow(((sky[:,:,z]-x[:,:,z])./sum(sky[:,:,z])),vmin = 0)
        e += 1
    end

    subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    cax = axes([0.85, 0.1, 0.025, 0.17])
    colorbar(cax=cax)
