

function [] = gen_quilting_map_for_proper(mp)

    %--PROPER initialization
    wl_dummy = 1e-6; %--dummy value needed to initialize wavelength in PROPER (meters)
    orderOfOps = 'XYZ';
    pixperact0 = 100; % Allowed values: 100 or 190
    inf_fn = sprintf('single_act_surface_in_m_bmc_2k_%dpixperact.fits', pixperact0);
    inf_sign = '+';
    cshift = 0;

    diamAtDM = 0.018511968; % beam diameter at DMs from the inside the DST2 PROPER model [meters]
    dx = diamAtDM / mp.P1.full.Nbeam; 
    bdf = 1;
    
    pixperact = mp.dm1.dm_spacing / dx;
    demagfac = pixperact0 / pixperact;
    
    % Find largest even integer component of demagfac to use for Fourier downsampling
    demagfacInt = floor(demagfac);
    if mod(demagfacInt, 2) == 1
        demagfacInt = demagfacInt - 1;
    end
    dxHD = dx / demagfacInt;
    Narray = ceil_even(mp.P1.full.Nbeam * 1.6);
    NarrayHD = demagfacInt * Narray;

    NactExtra = 20; % Have the quilting extend beyond the actuated area
    bm = prop_begin(NarrayHD*dxHD, wl_dummy, NarrayHD, bdf);

    dm = mp.dm1;
    [~, dm1SurfHD] = propcustom_dm_quilting(bm, ones(mp.dm1.Nact+NactExtra), dm.xc-cshift+NactExtra/2, dm.yc-cshift+NactExtra/2, ...
        dm.dm_spacing, 'XTILT', dm.xtilt, 'YTILT', dm.ytilt, 'ZTILT', dm.zrot, orderOfOps, ...
        'inf_sign', inf_sign, 'inf_fn', inf_fn);

    dm = mp.dm2;
    [~, dm2SurfHD] = propcustom_dm_quilting(bm, ones(mp.dm2.Nact+NactExtra), dm.xc-cshift+NactExtra/2, dm.yc-cshift+NactExtra/2, ...
        dm.dm_spacing, 'XTILT', dm.xtilt, 'YTILT', dm.ytilt, 'ZTILT', dm.zrot, orderOfOps, ...
        'inf_sign', inf_sign, 'inf_fn', inf_fn);

    % Fourier downsampling
    % Energy is conserved (except for the satellite spots that are cropped out.)
    dm1Surf =  real(ifftshift(ifft2(fftshift( pad_crop(fftshift(fft2(ifftshift(dm1SurfHD))), Narray))))) / demagfacInt;
    dm2Surf =  real(ifftshift(ifft2(fftshift( pad_crop(fftshift(fft2(ifftshift(dm2SurfHD))), Narray))))) / demagfacInt;
    
    if mp.flagPlot
        figure(651); imagesc(dm1SurfHD); axis xy equal tight; colorbar; title('DM1 Quilting'); drawnow;
        figure(652); imagesc(dm2SurfHD); axis xy equal tight; colorbar; title('DM2 Quilting'); drawnow;
        figure(661); imagesc(dm1Surf); axis xy equal tight; colorbar; title('DM1 Quilting'); drawnow;
        figure(662); imagesc(dm2Surf); axis xy equal tight; colorbar; title('DM2 Quilting'); drawnow;
    end
    %--Save out file
    import matlab.io.*

    fnDM1 = [mp.full.map_dir filesep 'dm1_surface_quilting.fits'];
    fnDM2 = [mp.full.map_dir filesep 'dm2_surface_quilting.fits'];

    fptrDM = fits.createFile(['!' fnDM1]); % '!' at front allows overwriting
    fits.createImg(fptrDM, 'double', size(dm1Surf));
    fits.writeImg(fptrDM, dm1Surf);
    fits.writeKey(fptrDM, 'PIXSIZE', dx, 'meters per pixel');
    fits.closeFile(fptrDM);

    fptrDM = fits.createFile(['!' fnDM2]); % '!' at front allows overwriting
    fits.createImg(fptrDM, 'double', size(dm2Surf));
    fits.writeImg(fptrDM, dm2Surf);
    fits.writeKey(fptrDM, 'PIXSIZE', dx, 'meters per pixel');
    fits.closeFile(fptrDM);

%     % Visually check that the fake satellite spots from aliasing are gone:
%     padFac = 8;
%     dm1Surf = dm1Surf - mean(dm1Surf(:));
%     window1D = tukeywindow(Narray, 0.4);
%     window2D = window1D * window1D.';
%     dm1SurfWindowed = dm1Surf .* window2D;
%     dmSurfPad = pad_crop(dm1SurfWindowed, padFac*mp.P1.full.Nbeam);
%     Efoc = fftshift(fft2(ifftshift(dmSurfPad)));
%     Ifoc = abs(Efoc).^2;
%     Ifoc = Ifoc/max(Ifoc(:));
%     figure(14); imagesc(log10(Ifoc), [-9, 0]); axis xy equal tight; colorbar;  drawnow;
%     figure(15); imagesc(dm1SurfWindowed); axis xy equal tight; colorbar; drawnow;
%     disp('');

end