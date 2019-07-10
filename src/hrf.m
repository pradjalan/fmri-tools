function ret = hrf(tr, nsamples)


    t = 1:1:tr*nsamples; % MEASUREMENTS
    h = gampdf(t,6, 1000) + -.5*gampdf(t,10, 1000); % HRF MODEL in ms resolution
    h = h/max(h); % SCALE HRF TO HAVE MAX AMPLITUDE OF 1

    i = 1:nsamples;
    ret = h(((i-1)*tr)+1);

end
