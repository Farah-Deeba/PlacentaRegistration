clear all
close all
clc

%% select which placenta to work on

placenta_id = '005'; %selecting an example placenta
placenta_str = sprintf('D:\\SWAVE2.0\\SWAVE2.0_DATA_ORGANIZED\\placenta%s',(placenta_id));
placenta_directory = dir(placenta_str);
saveMode = 0; % set it 0 to select the fiducials yourself
[ResultMRI,RA_MRI] = RGB2MRI_RegistrationSelectFiducials(placenta_id, saveMode);
       
%% function for RGB and MRI registration
function [axial_mri,axial_RA] = RGB2MRI_RegistrationSelectFiducials(placenta_id, saveMode)
   
    case_name = sprintf('%s',(placenta_id));
    var_string = sprintf('placenta%s',case_name);

    Placenta = load('placentaGit');
    placenta_jpeg = imread('Placenta005.jpg');
    
    file_name ='1.2.840.113619.2.408.4738430.15806453.21836.1528216298.7'; %coronal MRI
    file_orientation = dicominfo(file_name).ImageOrientationPatient; 
    MRI = dicomread(file_name);
    MRIinfo = dicominfo(file_name);
    
    fig1 = figure; mri = imagesc(MRI), colormap gray
   
    fig2 = figure; rgb = imagesc(placenta_jpeg)

    
%%

    
    if saveMode == 0
        prompt = {'Distance between capsules from physical ruler:'};
        dlgtitle = 'Distance between capsules (mm)';
        dims = [1 70];
        distRuler = inputdlg(prompt,dlgtitle,dims);
    else
        distRuler = Placenta.(var_string).distRuler;
    end

    if saveMode == 0
        h = images.roi.AssistedFreehand(mri);
        draw(h);
        X =  h.Position;

        h1= images.roi.AssistedFreehand(mri);
        draw(h1);
        X1 = h1.Position;
        X = [X;X1];
        h1= images.roi.AssistedFreehand(mri);
        draw(h1);
        X2 = h1.Position;
        X = [X;X2];

       
        rotAngle = 90; % select an angle so that RGB image is aligned to the coronal MRI
      

        placenta_jpeg_rot = imrotate(placenta_jpeg,(rotAngle));
        fig3 = figure; imagesc(placenta_jpeg_rot); colormap jet

        h= images.roi.AssistedFreehand;
        draw(h);
        Y = h.Position;
        h = images.roi.AssistedFreehand;
        draw(h);
        Y1 =  h.Position;
        Y = [Y; Y1];
        h =  images.roi.AssistedFreehand;
        draw(h);
        Y2 =  h.Position;
        Y = [Y; Y2];
        close (fig3)
    else
        Y = Placenta.(var_string).fiducial_pos_rgb;
        X = Placenta.(var_string).fiducial_pos_mri;
        rotAngle = Placenta.(var_string).rgb_rot_angle;
        placenta_jpeg_rot = imrotate(placenta_jpeg,str2double((rotAngle)));
    end
    
    img1 = im2double(placenta_jpeg_rot );
    refim = im2double(MRI);

    % Set the options
    opt.method='rigid'; % use rigid registration
    opt.viz=1;
    opt.scale =0; 
    opt.rot = 1 ;
    opt.corresp=1;
    opt.tol = 1e-7;
    % show every iteration
    % registering Y to X
    [Transform,~]=cpd_register(Y,X,opt);

    % Now lets apply the found transformation to the image
    % Create a dense grid of points
    [M,N]=size(refim);
    [x,y]=meshgrid(1:N,1:M); 
    grid=[x(:) y(:)];

    % Transform the grid according to the estimated transformation
    T=cpd_transform(grid, Transform);

    % Interpolate the image
    Tx=reshape(T(:,1),[M N]);
    Ty=reshape(T(:,2),[M N]);
    result(:,:,1)=interp2(img1(:,:,1),Tx,Ty);   % evaluate image Y at the new gird (through interpolation)
    result(:,:,2)=interp2(img1(:,:,2),Tx,Ty);   % evaluate image Y at the new gird (through interpolation)
    result(:,:,3)=interp2(img1(:,:,3),Tx,Ty);   % evaluate image Y at the new gird (through interpolation)


    % Show the result
    figure('Name','Registration Results'); subplot 311, 
    subplot 121,imshow(refim,[]); title('MRI');
    subplot 122,imshow(result,[]); title('Prime Snap');

  

    RA = imref2d(size(MRI),MRIinfo.PixelSpacing(2),MRIinfo.PixelSpacing(1));
   
    figure('Name','Images in World Coordinate');
    subplot 121, imshow(MRI,RA,'DisplayRange',[]), title('coronal')
    subplot 122, imshow(result,RA,'DisplayRange',[]), title('RGB')
     string1 = sprintf('distance(ruler) = %s cm', cell2mat(distRuler));
    text(60, 110, string1,'Color','w')
   
    if saveMode == 0
        prompt = {'what the markers indicate from left to right?'};
        dlgtitle = 'marker identification';
        dims = [1 70];
        markerID = cell2mat(inputdlg(prompt,dlgtitle,dims)); % L for lesion and N for normal; The black dot indicates lesion
    else
        markerID = Placenta.(var_string).markerID;
    end
    
    if saveMode ==0
    % mark the required points
        f1 = figure;
        imshow(rgb2hsv(result),RA,'DisplayRange',[])

        [x,y] = ginput(2+length(markerID)); % mark the four points sequentially: (1) the center of the left pill, (ii) center of right pill, (iii) left ROI1 and (iv)right ROI2
        hold on
        scatter(x,y,'filled','b')
        dist_capsule = pdist([x(1) y(1); [x(2) y(2)]]);
        dist_1 = pdist([x(1) y(1); [x(3) y(3)]]);
        if length(markerID)>1, dist_2 = pdist([x(1) y(1); [x(4) y(4)]]); end
        if length(markerID)>2, dist_3 = pdist([x(1) y(1); [x(5) y(5)]]); end
        close(f1)
        f2 = figure;
        imshow(result,RA,'DisplayRange',[])
        imdistline(gca,[x(1) x(2)],[y(1) y(2)]);
        string1 = sprintf('distance(ruler) = %s cm', cell2mat(distRuler));
        text(x(1)-10, y(1)+10, string1,'Color','w')
   
    end
        
 

    %
    axial_mri = dicomread('1.2.840.113619.2.408.4738430.15806453.21836.1528216298.448');
    axialMRI_info = dicominfo('1.2.840.113619.2.408.4738430.15806453.21836.1528216298.448');
    axial_RA = imref2d(size(axial_mri),axialMRI_info.PixelSpacing(2),axialMRI_info.PixelSpacing(1));
    
  
    if saveMode == 0
        % mark the centers of the oil capsules
        f1 = figure; imshow(axial_mri,axial_RA,'DisplayRange',[])
        [x2,y2] = ginput(2);
        dist_capsule_mri = pdist([x2(1) y2(1); [x2(2) y2(2)]]);
        v = [x2(2),y2(2)]-[x2(1),y2(1)];
        u = v./norm(v);
        dist_11 = [x2(1) y2(1)] + dist_1.*u;
        x2(3) = dist_11(1); y2(3) = dist_11(2);
        if length(markerID)>1; dist_22 = [x2(1) y2(1)] + dist_2.*u;
        x2(4) = dist_22(1); y2(4) = dist_22(2);end
        if length(markerID)>2; dist_33 = [x2(1) y2(1)] + dist_3.*u;
        x2(5) = dist_33(1); y2(5) = dist_33(2);end
        close (f1)
    else
        x2 = Placenta.(var_string).mriAxFiducialPos(:,1);
        y2 = Placenta.(var_string).mriAxFiducialPos(:,2);
    end
   
    
    overlay_design = zeros(size(axial_mri));
    [xIntrinsic, yIntrinsic] = worldToIntrinsic(axial_RA,x2,y2);
    xIntrinsic = round(xIntrinsic);  yIntrinsic = round(yIntrinsic);
    overlay_design(yIntrinsic(3)-20:yIntrinsic(3)+20,(xIntrinsic(3))) = 1;
    if length(markerID)>1; overlay_design(yIntrinsic(3)-20:yIntrinsic(3)+20,(xIntrinsic(4))) = 1;end
    if length(markerID)>2; overlay_design(yIntrinsic(3)-20:yIntrinsic(3)+20,(xIntrinsic(5))) = 1;end
%     figure, imshow(overlay_design,axial_RA,'DisplayRange',[])
    if strcmp(markerID(1),'L'); str1 = sprintf('Lesion'); else, str1 = sprintf('Normal');end
    if length(markerID)>1 && strcmp(markerID(2),'L' ); str2 = sprintf('Lesion'); else, str2 = sprintf('Normal');end
    if length(markerID)>2 && strcmp(markerID(3),'L' ); str3 = sprintf('Lesion'); else, str3 = sprintf('Normal');end

    overlay_design = insertText(overlay_design,[round(x2(3)),10],str1);
    if length(markerID)>1;overlay_design = insertText(overlay_design,[round(x2(4)),10],str2);end
    if length(markerID)>2;overlay_design = insertText(overlay_design,[round(x2(5)),10],str3); end
    
   
    ax_mri = figure('Name','Registered MRI');
    imshow(im2double(axial_mri)+overlay_design,axial_RA,'DisplayRange',[.4 .6]);
    if strcmp(markerID(1),'L'); str1 = sprintf('Lesion'); else, str1 = sprintf('Normal');end
    if length(markerID)>1 && strcmp(markerID(2),'L' ); str2 = sprintf('Lesion'); else, str2 = sprintf('Normal');end
    if length(markerID)>2 && strcmp(markerID(3),'L' ); str3 = sprintf('Lesion'); else, str3 = sprintf('Normal');end

    
    ax_mri = figure('Name','MRI (annotated with distance between capsules)');
    imshow(im2double(axial_mri),axial_RA,'DisplayRange',[.4 .6]);
    imdistline(gca,[x2(1) x2(2)],[y2(1) y2(2)]);
    
  
   
  
    
   end


