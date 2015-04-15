####################################################################
############################# Plot ADMM ############################
####################################################################

figure(1)
    clf()
    subplot(3,2,4)
    colorbar(imshow(x[:,:,1]))
    title("reconstruction")
    subplot(3,2,1)
    colorbar(imshow(sky[:,:,1]))
    title("sky")
    # subplot(3,1,3)
    # plot(snr[:])
    # title("snr ratio")
    # ylabel("dB")
    # xlabel("iterations")
    subplot(3,2,3)
    colorbar(imshow(mypsf[:,:,1]))
    title("psf")
    subplot(3,2,2)
    colorbar(imshow(mydata[:,:,1]))
    title("sky convolué")
