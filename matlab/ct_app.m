function ct_app()
% CT_APP  Minimal interactive interface for the CT simulation engine.
%   Every control maps to one field of the ct_simulate parameter struct.

fig = uifigure('Name','CT Simulator','Position',[100 100 900 560]);
g = uigridlayout(fig, [1 2]);
g.ColumnWidth = {260, '1x'};

% ---------- left: controls -------------------------------------------
p = uipanel(g,'Title','Acquisition');
c = uigridlayout(p,[14 1]);
c.RowHeight = repmat({'fit'},1,14);

uilabel(c,'Text','Tube voltage (kVp)');
kvpDD = uidropdown(c,'Items',{'50','80','120'},'Value','80');

uilabel(c,'Text','Geometry');
geoDD = uidropdown(c,'Items',{'parallel','fan'},'Value','parallel');

uilabel(c,'Text','Reconstruction');
recDD = uidropdown(c,'Items',{'fbp','sirt'},'Value','fbp');

uilabel(c,'Text','Filter');
filDD = uidropdown(c,'Items',{'Ram-Lak','Shepp-Logan','Cosine','Hamming','Hann'}, ...
                   'Value','Hann');

uilabel(c,'Text','Detail material');
matDD = uidropdown(c,'Items',{'none','water','pmma','al'},'Value','none');

statLbl = uilabel(c,'Text','Incident photons/pixel: 1e6');
statSld = uislider(c,'Limits',[3 8],'Value',6, ...
                   'MajorTicks',3:8,'MajorTickLabels',{'1e3','1e4','1e5','1e6','1e7','1e8'});
statSld.ValueChangedFcn = @(s,~) set(statLbl,'Text', ...
    sprintf('Incident photons/pixel: 1e%d', round(s.Value)));

runBtn = uibutton(c,'Text','Run simulation','BackgroundColor',[0.20 0.45 0.70], ...
                  'FontColor','w','FontWeight','bold');

% ---------- right: image + readout -----------------------------------
rp = uipanel(g,'Title','Reconstruction');
rg = uigridlayout(rp,[2 1]);
rg.RowHeight = {'1x', 'fit'};
ax = uiaxes(rg);
readout = uilabel(rg,'Text','Press Run to simulate.','FontSize',13);

runBtn.ButtonPushedFcn = @(~,~) doRun();

% ---------- callback ---------------------------------------------------
    function doRun()
        runBtn.Enable = 'off';
        readout.Text = 'Running...';  drawnow;
        try
            q = struct();
            q.kvp          = str2double(kvpDD.Value);
            q.geometry     = geoDD.Value;
            q.recon_method = recDD.Value;
            q.recon_filter = filDD.Value;
            q.n0_per_pixel = 10^round(statSld.Value);
            if strcmp(matDD.Value,'none')
                q.det_material = '';
            else
                q.det_material = matDD.Value;
            end
            if strcmp(q.recon_method,'sirt'), q.n_angles = 60; end

            r = ct_simulate(q);

            imagesc(ax, r.recon); axis(ax,'image','off'); colormap(ax,'gray');
            colorbar(ax);
            title(ax, sprintf('%s, %d kVp, %s', q.geometry, q.kvp, q.recon_method));

            txt = sprintf(['\\mu centre = %.4f    \\mu edge = %.4f    ' ...
                           'cupping = %.2f%%    noise sd = %.4f'], ...
                          r.mu_centre, r.mu_edge, r.cupping, r.noise_sd);
            if ~isnan(r.contrast)
                txt = [txt sprintf('    contrast = %.1f%%', r.contrast)];
            end
            readout.Text = txt;
            readout.Interpreter = 'tex';
        catch ME
            readout.Text = ['Error: ' ME.message];
        end
        runBtn.Enable = 'on';
    end
end