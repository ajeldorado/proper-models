

function [] = gen_quilting_map_for_proper(mp)

    %--PROPER initialization
    wl_dummy = 1e-6; %--dummy value needed to initialize wavelength in PROPER (meters)
    orderOfOps = 'XYZ';
    inf_fn = 'single_act_surface_in_m_bmc_2k_100pixperact.fits';
    inf_sign = '+';
    cshift = 0;

    % mp.full.gridsize = 1024;
    % mp.P1.full.Nbeam = 350;

    diamAtDM = 0.018511968; % beam diameter at DMs from the PROPER inside the model [meters]
    dx = diamAtDM / mp.P1.full.Nbeam; 
    Narray = ceil_even(2*mp.P1.full.Nbeam);
    bdf = 1;%mp.P1.full.Nbeam / Narray;

    NactExtra = 20; % Have the quilting extend beyond the actuated area
    bm = prop_begin(Narray*dx, wl_dummy, Narray, bdf);

    dm = mp.dm1;
    [~, dm1Surf] = propcustom_dm_quilting(bm, ones(mp.dm1.Nact+NactExtra), dm.xc-cshift+NactExtra/2, dm.yc-cshift+NactExtra/2, ...
        dm.dm_spacing, 'XTILT', dm.xtilt, 'YTILT', dm.ytilt, 'ZTILT', dm.zrot, orderOfOps, ...
        'inf_sign', inf_sign, 'inf_fn', inf_fn);

    dm = mp.dm2;
    [~, dm2Surf] = propcustom_dm_quilting(bm, ones(mp.dm2.Nact+NactExtra), dm.xc-cshift+NactExtra/2, dm.yc-cshift+NactExtra/2, ...
        dm.dm_spacing, 'XTILT', dm.xtilt, 'YTILT', dm.ytilt, 'ZTILT', dm.zrot, orderOfOps, ...
        'inf_sign', inf_sign, 'inf_fn', inf_fn);

    if mp.flagPlot
        figure(651); imagesc(dm1Surf); axis xy equal tight; colorbar; title('DM1 Quilting'); drawnow;
        figure(652); imagesc(dm2Surf); axis xy equal tight; colorbar; title('DM2 Quilting'); drawnow;
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

    % padFac = 8;
    % dmSurfPad = pad_crop(dm1Surf, padFac*mp.P1.full.Nbeam);
    % Efoc = fftshift(fft2(ifftshift(dmSurfPad)));
    % Ifoc = abs(Efoc).^2;
    % Ifoc = Ifoc/max(Ifoc(:));
    % figure(14); imagesc(log10(Ifoc), [-5, 0]); axis xy equal tight; colorbar;


end