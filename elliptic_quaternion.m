classdef elliptic_quaternion

    %Identifies the elliptic quaternion
    properties
        r_part
        i_part
        j_part
        k_part
        p
    end

    methods (Access=public)

        %Constructor
        function obj = elliptic_quaternion(r_part,i_part,j_part,k_part,p)
            if ~isequal(size(r_part), size(i_part), size(j_part), size(k_part))|| p>=0
                error('ELLIPTIC:constructor','r_part,i_part,j_part and k_part must be same size and p must be negative')
            else
                obj.r_part = r_part;
                obj.i_part = i_part;
                obj.j_part = j_part;
                obj.k_part=k_part;
                obj.p=p;
            end
        end

        %Returns the real part of the elliptic quaternion
        function r = elliptic_q_getrealpart(obj)
            r = obj.r_part;
        end

        %Returns the imaginary part of the elliptic quaternion
        function [i_part, j_part, k_part] = elliptic_q_getimagparts(obj)
            i_part=obj.i_part;
            j_part=obj.j_part;
            k_part=obj.k_part;
        end

        %Returns the elliptic quaternion (or elliptic quaternion matrix) given
        % in ijk -form in $e_1$-$e_2$ form.
        function [e1_part,e2_part]=e1e2form(obj)
            e1_part=elliptic_number(obj.r_part+obj.j_part,obj.i_part+obj.k_part,obj.p);
            e2_part=elliptic_number(obj.r_part-obj.j_part,obj.i_part-obj.k_part,obj.p);
        end

        %Returns the elliptic quaternion (or elliptic quaternion matrix) given
        % in $e_1$-$e_2$ form in ijk form.
        function  obj=ijkform(e1_part,e2_part)
            rpart=(e1_part.r_part+e2_part.r_part)/2;
            jpart=(e1_part.r_part-e2_part.r_part)/2;
            ipart=(e1_part.i_part+e2_part.i_part)/2;
            kpart=(e1_part.i_part-e2_part.i_part)/2;
            obj=elliptic_quaternion(rpart,ipart,jpart,kpart,e1_part.p);
        end

        %Returns the 1st conjugate of the elliptic quaternion (or elliptic quaternion matrix).
        function obj=elliptic_q_conjugate_1(obj1)
            obj=elliptic_quaternion(obj1.r_part,-obj1.i_part,obj1.j_part,-obj1.k_part,obj1.p);
        end
        %Returns the 2nd conjugate of the elliptic quaternion (or elliptic quaternion matrix).
        function obj=elliptic_q_conjugate_2(obj1)
            obj=elliptic_quaternion(obj1.r_part,obj1.i_part,-obj1.j_part,-obj1.k_part,obj1.p);
        end
        %Returns the 3rd conjugate of the elliptic quaternion (or elliptic quaternion matrix).
        function obj=elliptic_q_conjugate_3(obj1)
            obj=elliptic_quaternion(obj1.r_part,-obj1.i_part,-obj1.j_part,obj1.k_part,obj1.p);
        end

        %Returns the elliptic quaternion (or elliptic quaternion matrix) times -1.
        function obj = elliptic_q_minus(obj1)
            obj = elliptic_quaternion(-obj1.r_part, -obj1.i_part, -obj1.j_part,-obj1.k_part,obj1.p);
        end

        %Returns the sum of two elliptic quaternions (or elliptic quaternion matrices).
        function obj=elliptic_q_addition(obj1,obj2)
            obj=elliptic_quaternion(obj1.r_part+obj2.r_part,obj1.i_part+obj2.i_part,obj1.j_part+obj2.j_part,obj1.k_part+obj2.k_part,obj1.p);
        end
        %Returns the difference of two elliptic quaternions (or elliptic quaternion matrices).
        function obj=elliptic_q_subtraction(obj1,obj2)
            obj=elliptic_q_addition(obj1,elliptic_q_minus(obj2));
        end

        %Returns the product of two elliptic quaternions (or elliptic quaternion matrices).
        function obj=elliptic_q_product(obj1,obj2)
            [a_e1_part,a_e2_part]=e1e2form(obj1);
            [b_e1_part,b_e2_part]=e1e2form(obj2);
            obj=ijkform(elliptic_n_product(a_e1_part,b_e1_part),elliptic_n_product(a_e2_part,b_e2_part));
        end

        %Returns the lambda scalar multiple of the elliptic quaternion (or elliptic quaternion matrix).
        function obj=elliptic_q_scalar(lamda,obj1)
            obj=elliptic_quaternion(lamda*obj1.r_part,lamda*obj1.i_part,lamda*obj1.j_part,lamda*obj1.k_part,obj1.p);
        end

        % Checks if the elliptic quaternion is zero. If it is zero,
        % it returns 1, otherwise it returns 0.
        function [output]=elliptic_q_isazerodivisor(obj)
            [e1_part,e2_part]=e1e2form(obj);
            if elliptic_n_isazero(e1_part)==1 || elliptic_n_isazero(e2_part)==1
                output=true;
            else
                output=false;
            end
        end

        %Returns inverse of the elliptic quaternion (or elliptic quaternion matrix).
        function obj=elliptic_q_inverse(obj1)
            [e1_part,e2_part]=e1e2form(obj1);
            e1_part=elliptic_n_inverse(e1_part);
            e2_part=elliptic_n_inverse(e2_part);
            obj=ijkform(e1_part,e2_part);
        end

        %Returns the division of two elliptic quaternion (or elliptic quaternion matrix).
        function obj = elliptic_q_divide(obj1,obj2)
            if elliptic_q_isazerodivisor(obj2)==1
                fprintf('Error executing code: Division by zero');
            else
                obj = elliptic_q_product(obj1,elliptic_q_inverse(obj2));
            end
        end

        %Checks the equality of two elliptic quaternions (or elliptic quaternion matrices).
        function [output]=elliptic_q_equal(obj1,obj2)
            if (obj1.r_part==obj2.r_part && obj1.i_part==obj2.i_part && obj1.j_part==obj2.j_part && obj1.k_part==obj2.k_part && obj1.p==obj2.p)
                output=true;
            else
                output=false;
            end
        end

        %Returns the power of elliptic quaternion (or elliptic quaternion matrix).
        function obj=elliptic_q_power(obj1,n)
            [e1_part,e2_part]=e1e2form(obj1);
            e1_part=elliptic_n_power(e1_part,n);
            e2_part=elliptic_n_power(e2_part,n);
            obj=ijkform(e1_part,e2_part);
        end

        %Returns the norm of the elliptic quaternion (or elliptic quaternion matrix).
        function q_norm=elliptic_q_norm(obj)
            q_norm=sqrt((norm(obj.r_part,'fro'))^2-obj.p*(norm(obj.i_part,'fro'))^2+(norm(obj.j_part,'fro'))^2-obj.p*(norm(obj.k_part,'fro'))^2);
        end

        %Returns the transpose of the elliptic quaternion matrix.
        function obj=elliptic_q_transpose(obj1)
            [e1_part,e2_part]=e1e2form(obj1);
            e1_part=elliptic_n_transpose(e1_part);
            e2_part=elliptic_n_transpose(e2_part);
            obj=ijkform(e1_part,e2_part);
        end

        %Returns the conjugate transpose of the elliptic quaternion matrix.
        function obj=elliptic_q_hermitianconjugate(obj1)
            [e1_part,e2_part]=e1e2form(obj1);
            e1_part=elliptic_n_hermitianconjugate(e1_part);
            e2_part=elliptic_n_hermitianconjugate(e2_part);
            obj=ijkform(e1_part,e2_part);
        end

        %Returns the ith row elements of the elliptic quaternion matrix
        function obj=elliptic_q_i_th_eqrow(obj1,i)
            [e1_part,e2_part]=e1e2form(obj1);
            obj=ijkform(elliptic_n_i_th_erow(e1_part,i),elliptic_n_i_th_erow(e2_part,i));
        end

        %Returns the ith column elements of the elliptic quaternion matrix
        function obj=elliptic_q_i_th_eqcolumn(obj1,i)
            [e1_part,e2_part]=e1e2form(obj1);
            obj=ijkform(elliptic_n_i_th_ecolumn(e1_part,i),elliptic_n_i_th_ecolumn(e2_part,i));
        end

        %Returns the ith row and jth column element of the elliptic complex matrix
        function obj=elliptic_q_index(obj1,i,j)
            [e1_part,e2_part]=e1e2form(obj1);
            obj=ijkform(elliptic_n_index(e1_part,i,j),elliptic_n_index(e2_part,i,j));
        end

        %Indexes the elliptic quaternion matrix.
        function obj=elliptic_q_indexing(obj1,r1,r2,c1,c2)
            [e1_part,e2_part]=e1e2form(obj1);
            obj1=elliptic_n_index(e1_part,r1,r2,c1,c2);
            obj2=elliptic_n_index(e2_part,r1,r2,c1,c2);
            obj=ijkform(obj1,obj2);
        end

        %Returns diagonal matrix obj1 of  eigenvalues and  matrix obj2
        % whose columns are the corresponding eigenvectors, so that
        % $obj*obj2 = obj2*obj1
        function [obj1,obj2]=elliptic_q_eig(obj)
            [e1_part,e2_part]=e1e2form(obj);

            eta_A1=ellipticcomplex2complexnumber(e1_part);
            [EigVec1,EigVal1]=eig(eta_A1);
            obj_1=complexnumber2ellipticcomplex(e1_part,EigVec1);
            obj_2=complexnumber2ellipticcomplex(e1_part,EigVal1);

            eta_A2=ellipticcomplex2complexnumber(e2_part);
            [EigVec2,EigVal2]=eig(eta_A2);
            obj_3=complexnumber2ellipticcomplex(e1_part,EigVec2);
            obj_4=complexnumber2ellipticcomplex(e1_part,EigVal2);

            obj1=ijkform(obj_1,obj_3);
            obj2=ijkform(obj_2,obj_4);
        end

        %Returns the pseudo inverse of the elliptic quaternion matrix.
        function obj=elliptic_q_pseudoinverse(obj1)
            [e1_part,e2_part]=e1e2form(obj1);
            obj=ijkform(elliptic_n_pinv(e1_part),elliptic_n_pinv(e2_part));
        end

        %[obj1,obj2,obj3] = elliptic_q_svd(obj) performs a singular value
        % decomposition of elliptic quaternion matrix obj,
        % such that obj = obj1*obj2*elliptic_q_hermitianconjugate(obj3).
        function [obj1,obj2,obj3]=elliptic_q_svd(obj)
            [e1_part,e2_part]=e1e2form(obj);
            [U1,S1,V1]=elliptic_n_svd(e1_part);
            [U2,S2,V2]=elliptic_n_svd(e2_part);
            obj1=ijkform(U1,U2);
            obj2=ijkform(S1,S2);
            obj3=ijkform(V1,V2);
        end

        %Returns the rank of the elliptic quaternion matrix
        function q_rank=elliptic_q_rank(obj)
            [e1_part,e2_part]=e1e2form(obj);
            q_rank=max(elliptic_n_rank(e1_part),elliptic_n_rank(e2_part));
        end

        %Returns the least squares solution of the elliptic  quaternion matrix equation obj1*obj=obj2
        % and the least square error.
        function [obj, e_error]=elliptic_q_lss(obj1,obj2)
            obj=elliptic_q_product(elliptic_q_pseudoinverse(obj1),obj2);
            e_error=elliptic_q_minus(elliptic_q_minus(elliptic_q_product(obj1,obj),obj2));
        end

        %Returns the optimal p value for the least squares solution of the
        %elliptic  quaternion matrix equation obj_A*obj=obj_b.
        function optimal_p_for_lss(obj_A,obj_b,start_p,end_p,step_size)
            x_axes=start_p:step_size:end_p;
            y_axes=zeros(length(x_axes));
            for i=1:length(x_axes)
                obj_A.p=x_axes(i);
                obj_b.p=x_axes(i);
                [~, e_error]=elliptic_q_lss(obj_A,obj_b);
                y_axes(i)=e_error;
            end
            [min_error, min_index] = min(y_axes);
            optimal_p = x_axes(min_index);
            figure
            plot(optimal_p, min_error, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
            hold on
            plt=plot(x_axes, y_axes);
            grid on;
            grid(gca, 'minor');
            set(gca, 'GridLineStyle', '-');
            set(gca, 'GridColor', [0.5, 0.5, 0.5]);
            plt.Color='b';
            plt.Marker='o';
            plt.LineStyle='-';
            plt.LineWidth=1.5;
            hold off
            xlabel('p');
            ylabel('error');
            optimal_p_value="optimal p value:"+num2str(optimal_p);
            legend(optimal_p_value, 'Location', 'Best');
        end

        %Returns the elliptic quaternion matrix representation of the image.
        function obj=image2elliptic_q(obj1)
            %close all
            [filename, pathname] = uigetfile({'*.png';'*.jpeg';'*.jpg'},'Select MRI');
            inputimage=strcat(pathname, filename);
            inputimage = imread(inputimage);
            inputimage =im2double(inputimage);
            % figure
            %imshow(inputimage)
            % title('Original Image');
            R=inputimage(:,:,1);
            G=inputimage(:,:,2);
            B=inputimage(:,:,3);
            obj=elliptic_quaternion(zeros(size(R)),R,G,B,obj1.p);
        end

        %Returns the elliptic quaternion matrix representation of the image.
        function obj=image2elliptic_q_1(obj1,I)
            %close all

            inputimage =im2double(I);
            % figure
            %imshow(inputimage)
            % title('Original Image');
            R=inputimage(:,:,1);
            G=inputimage(:,:,2);
            B=inputimage(:,:,3);
            obj=elliptic_quaternion(zeros(size(R)),R,G,B,obj1.p);
        end

        %Returns the RGB space representation of an image represented by an elliptic quaternion matrix.
        function I=elliptic_q2image(obj)
            R=obj.i_part;
            G=obj.j_part;
            B=obj.k_part;
            I=cat(3,R,G,B);
            I=im2double(I);
        end

        %Returns the desired eigenimage of the image
        %obj is the image matrix, and h is the eigenimage number
        function  image=elliptic_q_eigenimage(obj,h)
            [obj1,~,obj3]=elliptic_q_svd(obj);
            c=zeros(size(obj.r_part,1),size(obj.r_part,2));
            image=elliptic_quaternion(c,c,c,c,obj.p);
            eimage=elliptic_q_product(i_th_eqcolumn(obj1,h),elliptic_q_hermitianconjugate(i_th_eqcolumn(obj3,h)));
            newimage=elliptic_to_image(eimage);
            absolute_image = abs(newimage);
            min_value = min(absolute_image(:));
            max_value = max(absolute_image(:));
            normalized_image = (absolute_image - min_value) / (max_value - min_value);
            imshow(normalized_image);
            title(sprintf('%d st eigenimage for airplane.png', h));
        end

        % It shows the change of singular values of the image in the elliptical quaternion matrix space.
        function  obj1=elliptic_q_singularvalues(obj)
            close all
            [~,obj1,~]=elliptic_q_svd(obj);
            y_real=zeros(elliptic_q_rank(obj));
            y_jpart=zeros(elliptic_q_rank(obj));
            x=zeros(1,elliptic_q_rank(obj));
            for j=1:elliptic_q_rank(obj)
                y_jpart(j)=log10(obj1.j_part(j,j));
                y_real(j)=log10(obj1.r_part(j,j));
                x(j)=j;
            end
            subplot(1,2,2);
            plot(x,y_jpart,'lineWidth',2);
            title('j-part of singular values of the barbara.png')
            ylabel('log10 scale')
            subplot(1,2,1);
            plot(x,y_real,'r','lineWidth',2);
            title('r-part of singular values of the barbara.png')
            ylabel('log10 scale')
        end

        %It compresses the image according to the desired k value in the elliptical quaternion matrix space.
        %k (Type: integer): The number of singular values to retain for the reconstruction.
        function [I_c,I_o,E_c,E_o]=elliptic_q_image_reconstruction(obj,k)
            E_o=obj;
            I_o=elliptic_q2image(obj);
            [obj1,obj2,obj3]=elliptic_q_svd(obj);
            obj1=elliptic_q_m_indexing(obj1,1,size(obj.r_part,1),1,k);
            obj2=elliptic_q_m_indexing(obj2,1,k,1,k);
            obj3=elliptic_q_m_indexing(obj3,1,size(obj.r_part,2),1,k);
            image=elliptic_q_product(elliptic_q_product(obj1,obj2),elliptic_q_hermitianconjugate(obj3));
            E_c=image;
            I_c=elliptic_to_image(image);
            figure
            imshow(I_c);
            %PSNR=psnr(I_c,I_o);
            %mse=immse(I_c,I_o);
            % figure
            % subplot(1,2,1);
            % imshow(I_o);
            % title('Original image');
            % subplot(1,2,2)
            % imshow(I_c);
            % title('Compressed image for k=50 and p=-0.5')
            %fprintf('%s%1.5f','PSNR:',PSNR);
            %fprintf('%s%1.5f','MSE',immse(I_c,I_o));
            %figure;
            %imshow(I_c);
            % title('compressed image');
            % figure;
            % imshow(I_o);
            % title('original image');
        end

        %Returns the optimal p (for psnr) value for image compression performed in elliptic quaternion matrix space.
        function optimal_p_value_for_psnr=elliptic_q_optimal_p_for_psnr(obj)
            %close all
            %t=timeit;
            x_axes=10:20:150; %k_value
            y_axes=-5:0.1:-0.1; %p_value
            [X,Y]=meshgrid(x_axes,y_axes);
            Z=zeros(length(y_axes),length(x_axes));
            for i=1:length(y_axes)
                for j=1:length(x_axes)
                    obj.p=y_axes(i);
                    [I_c,I_o,~,~]=elliptic_q_image_reconstruction(obj,x_axes(j));
                    Z(i,j)= psnr(I_c,I_o);
                end
            end
            figure;
            subplot(1,2,1)
            colormap('jet');
            surface(X, Y, Z, 'FaceAlpha', 0.5);
            grid on;
            grid(gca, 'minor');
            set(gca, 'GridLineStyle', '-');
            set(gca, 'GridColor', [0.5, 0.5, 0.5]);
            %colorbar;
            view(3)
            xlabel('K value');
            ylabel('p value');
            zlabel('psnr value');
            [max_values, max_indices] = max(Z, [], 1);
            hold on;
            scatter3(X(1, :), Y(max_indices), max_values, 50, 'r', 'filled');
            title('(a)')
            hold off;
            %legend(' ','Optimal p values');
            optimal_p_value_for_psnr=y_axes(max_indices);
        end

        % Returns the optimal p (for mse) value for image compression performed in elliptic quaternion matrix space.
        function optimal_p_value_for_mse=elliptic_q_optimal_p_for_mse(obj)
            %close all
            %t=timeit;
            x_axes=10:20:150; %k_value
            y_axes=-5:0.1:-0.1; %p_value
            [X,Y]=meshgrid(x_axes,y_axes);
            Z=zeros(length(y_axes),length(x_axes));
            for i=1:length(y_axes)
                for j=1:length(x_axes)
                    obj.p=y_axes(i);
                    [I_c,I_o,~,~]=elliptic_q_image_reconstruction(obj,x_axes(j));
                    Z(i,j)= immse(I_c,I_o);
                end
            end
            subplot(1,2,2)
            colormap(parula(5));
            surface(X, Y, Z, 'FaceAlpha', 0.5);
            grid on;
            grid(gca, 'minor');
            set(gca, 'GridLineStyle', '-');
            set(gca, 'GridColor', [0.5, 0.5, 0.5]);
            %colorbar;
            view(3)
            xlabel('K value');
            ylabel('p value');
            zlabel('mse value');
            [min_values, min_indices] = min(Z, [], 1);
            hold on;
            scatter3(X(1, :), Y(min_indices), min_values, 50, 'r', 'filled');
            title('(b)')
            hold off;
            %legend(' ','Optimal p values');
            optimal_p_value_for_mse=y_axes(min_indices);
        end

        %Original host image obj1
        %Original watermark obj2
        %Watermarked host image obj
        function obj=eliptic_imresize(obj1,row,col)
            obj.r_part=imresize(obj1.r_part,[row col]);
            obj.i_part=imresize(obj1.i_part,[row col]);
            obj.j_part=imresize(obj1.j_part,[row col]);
            obj.k_part=imresize(obj1.k_part,[row col]);
            obj=elliptic_quaternion(obj1.r_part,obj1.i_part,obj1.j_part,obj1.k_part,obj.p);
        end

        %The embedding_elliptic_watermarking function embeds a watermark into an object using elliptic quaternion matrix algebra.
        %obj1: The primary object to be watermarked.
        %obj2: The watermark object to be embedded.

        function [water_marked_image,U_w,V_w]=elliptic_q_embedding_watermarking(obj1,obj2,alpha)
            M=size(obj1.r_part,1);
            N=size(obj1.r_part,2);
            obj2=eliptic_imresize(obj2,M,N);
            [U,S,V]=elliptic_q_svd(obj1);
            [U_w,S_w,V_w]=elliptic_q_svd(addition(S,scalar(alpha,obj2)));
            water_marked_image=elliptic_q_product(elliptic_q_product(U,S_w),elliptic_q_hermitianconjugate(V));
            % I=elliptic_to_image(water_marked_image);
            % w_psnr=psnr(elliptic_q2image(original_host_image),elliptic_to_image(water_marked_image))
            % w_mse=immse(elliptic_q2image(original_host_image),elliptic_to_image(water_marked_image))
            % figure
            % subplot(1,3,1);
            % imshow(I)
            % title('Original Image')
            % subplot(1,3,2);
            % imshow(elliptic_to_image(original_water_mark))
            % title('Watermark')
            % subplot(1,3,3);
            % imshow(elliptic_to_image(water_marked_image))
            % title('Watermarked Image')
        end

        %The extraction_elliptic_watermarking  function extracts a watermark from a watermarked object using the elliptic quaternion matrix algebra.
        %obj1: The watermarked object from which the watermark is to be extracted. \\ &
        %obj2: The original, unwatermarked object.
        function water_mark=elliptic_q_extraction_watermarking(obj1,obj2,U_w,V_w,alpha)
            [~,S_w,~]=elliptic_q_svd(obj1);
            C=elliptic_q_product(elliptic_q_product(U_w,S_w),elliptic_q_hermitianconjugate(V_w));
            [~,S,~]=elliptic_q_svd(obj2);
            water_mark=scalar(1/alpha,elliptic_q_minus(C,S));
            I=elliptic_to_image(water_mark);
            figure
            imshow(I)
        end

        %Finds the optimal p (for psnr) for watermarking using the elliptic quaternion matrix algebra
        %obj1 (Type: Object): The primary object to be watermarked.
        %obj2 (Type: Object):   The watermark object to be embedded.
        function elliptic_q_watermarking_optimal_p_for_psnr(obj1,obj2,alpha)
            x_axes=-0.2:0.001:-0.001;
            y_axes=zeros(1,length(x_axes));
            for i=1:length(x_axes)
                obj1.p=x_axes(i);
                obj2.p=x_axes(i);
                [water_marked,~,~]=elliptic_q_embedding_watermarking(obj1,obj2,alpha);
                water_marked_image=elliptic_to_image(water_marked);
                original_host_image_I=elliptic_q2image(obj1);
                y_axes(i)=psnr(original_host_image_I,water_marked_image);
                close all
            end
            [max_psnr, max_index] = max(y_axes);
            optimal_p = x_axes(max_index);
            figure
            plot(optimal_p, max_psnr, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
            hold on
            plot(x_axes, y_axes,"b");
            grid on;
            grid(gca, 'minor');
            set(gca, 'GridLineStyle', '-');
            set(gca, 'GridColor', [0.5, 0.5, 0.5]);
            hold off
            % title('mse vs p');
            xlabel('p-values');
            ylabel('PSNR');
            optimal_p_value="optimal p value:"+num2str(optimal_p);
            legend(optimal_p_value, 'Location', 'Best');

        end

        %Finds the optimal p (for psnr) for watermarking using the elliptic quaternion matrix algebra
        %obj1 (Type: Object): The primary object to be watermarked.
        %obj2 (Type: Object):   The watermark object to be embedded.
        function elliptic_q_watermarking_optimal_p_for_mse(obj1,obj2,alpha)
            x_axes=-1:0.001:-0.001;
            y_axes=zeros(1,length(x_axes));
            for i=1:length(x_axes)
                obj1.p=x_axes(i);
                obj2.p=x_axes(i);
                [water_marked,~,~]=elliptic_q_embedding_watermarking(obj1,obj2,alpha);

                water_marked_image=elliptic_to_image(water_marked);
                original_host_image_I=elliptic_q2image(obj1);
                y_axes(i)=immse(original_host_image_I,water_marked_image);
                close all
            end
            [min_psnr, min_index] = min(y_axes);
            optimal_p = x_axes(min_index);
            figure
            plot(optimal_p, min_psnr, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
            hold on
            plot(x_axes, y_axes,"b");
            grid on;
            grid(gca, 'minor');
            set(gca, 'GridLineStyle', '-');
            set(gca, 'GridColor', [0.5, 0.5, 0.5]);
            hold off
            % title('mse vs p');
            xlabel('p');
            ylabel('MSE');
            optimal_p_value="optimal p value:"+num2str(optimal_p);
            legend(optimal_p_value, 'Location', 'Best');

        end

    end %methods

end %classdef