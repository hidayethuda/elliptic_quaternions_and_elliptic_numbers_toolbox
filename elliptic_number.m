classdef elliptic_number < elliptic_quaternion

    %Static Methods p-trigonometric functions
    methods (Static)
        function value=sinp(tetha,p)
            value=(1/sqrt(abs(p)))*sin(tetha*sqrt(abs(p)));
        end

        function value=cosp(tetha,p)
            value=cos(tetha*sqrt(abs(p)));
        end

        function value=tanp(tetha,p)
            value=sinp(tetha,p)/cosp(tetha,p);
        end

    end %static methods

    %Methods of the class of elliptic numbers
    methods (Access=public)

        % Constructor
        function obj = elliptic_number(r_part,i_part,p)
            obj@elliptic_quaternion(r_part,i_part,zeros(size((r_part))),zeros((size(r_part))),p);
        end

        %Returns the complex number (or matrix) to which the 
        %elliptic number  (or matrix) is isomorphic.
        function eta=ellipticcomplex2complexnumber(obj)
            eta=obj.r_part+1i*sqrt(abs(obj.p))*obj.i_part;
        end

        %Returns the elliptic number (or matrix) to which the complex 
        %number (or matrix) is isomorphic.
        function obj=complexnumber2ellipticcomplex(obj1,complex)
            obj = elliptic_number(real(complex),(1/sqrt(abs(obj1.p)))*imag(complex),obj1.p);
        end

        %Returns the sum of two elliptic numbers (or matrices).
        function obj=elliptic_n_addition(obj1,obj2)
            obj=elliptic_number(obj1.r_part+obj2.r_part,obj1.i_part+obj2.i_part,obj1.p);
        end

        %Returns the difference of two elliptic numbers (or matrices).
        function obj=elliptic_n_subtraction(obj1,obj2)
            obj=elliptic_number(obj1.r_part-obj2.r_part,obj1.i_part-obj2.i_part,obj1.p);
        end

        %Returns the product of two elliptic numbers (or matrices).
        function obj=elliptic_n_product(obj1,obj2)
            obj=elliptic_number(obj1.r_part*obj2.r_part+obj1.p*obj1.i_part*obj2.i_part,obj1.r_part*obj2.i_part+obj1.i_part*obj2.r_part,obj1.p);
        end

        %Returns inverse of the elliptic numbers (or matrices).
        function obj = elliptic_n_inverse(obj1)
            eta_a=ellipticcomplex2complexnumber(obj1);
            eta_a_inv=inv(eta_a);
            obj=complexnumber2ellipticcomplex(obj1,eta_a_inv);
        end

        % Checks if the elliptic number is zero. If it is zero,
        % it returns 1, otherwise it returns 0.
        function output=elliptic_n_isazero(obj)
            if (obj.r_part==0 && obj.i_part==0)
                output=true;
            else
                output=false;
            end
        end
    
        %Returns the division of two elliptic numbers (or matrices).
        function obj = elliptic_n_divide(obj1,obj2)
            if elliptic_n_isazero(obj2)==1
                fprintf('Error executing code: Division by zero');
            else
                obj = elliptic_n_product(obj1,inverse(obj2));
            end
        end

        %Returns distance between two elliptic numbers.
        function distance=elliptic_n_distance(obj1,obj2)
            distance=sqrt((obj1.r_part-obj2.r_part).^2-obj1.p*(obj1.i_part-obj2.i_part).^2);
        end

        %Returns argument of the elliptic numbers.
        function argument = elliptic_n_argument(obj)
            if elliptic_n_isazero(obj)==0
                argument=(1/sqrt(abs(obj.p)))*atan((obj.i_part/obj.r_part)*sqrt(abs(obj.p)));
            end
        end

        %Returns the power of elliptic numbers.
        function obj=elliptic_n_power(obj1,n)
            obj=elliptic_number(norm(obj1).^n*obj1.cosp(n*elliptic_n_argument(obj1),obj1.p),norm(obj1).^n*obj1.sinp(n*elliptic_n_argument(obj1),obj1.p),obj1.p);
        end

        %Indexes the elliptic matrices. Stated in other words, gives a chance 
        %to choose which part of the elliptic complex matrix to use.
        function obj=elliptic_n_indexing(obj1,r1,r2,c1,c2)
            obj=elliptic_number(obj1.r_part(r1:r2,c1:c2),obj1.i_part(r1:r2,c1:c2),obj1.p);
        end

        %Returns the transpose of the elliptic matrices
        function obj=elliptic_n_transpose(obj1)
            eta_a=ellipticcomplex2complexnumber(obj1);
            eta_a_transpose=transpose(eta_a);
            obj=complexnumber2ellipticcomplex(obj1,eta_a_transpose);
        end

        %Returns the conjugate of the elliptic matrices
        function obj=elliptic_n_conjugate(obj1)
            eta_a=ellipticcomplex2complexnumber(obj1);
            eta_a_conjugate=conj(eta_a);
            obj=complexnumber2ellipticcomplex(obj1,eta_a_conjugate);
        end

        %Returns the Hermitian conjugate of the elliptic matrices
        function obj=elliptic_n_hermitianconjugate(obj1)
            obj=elliptic_transpose(elliptic_n_conjugate(obj1));
        end

        %Returns eigenvalues and eigenvectors of the elliptic matrices
        %obj1: eigenvectors, obj2: eigenvalues
        function [obj1,obj2]=elliptic_n_eig(obj)
            eta_A=ellipticcomplex2complexnumber(obj);
            [EigVec,EigVal]=eig(eta_A);
            obj1=complexnumber2ellipticcomplex(obj,EigVec);
            obj2=complexnumber2ellipticcomplex(obj,EigVal);
        end

        %Finds SVD parameters of the elliptic matrices
        %obj1:U , obj2:S, obj3: 
        function [obj1,obj2,obj3]=elliptic_n_svd(obj)
            eta=ellipticcomplex2complexnumber(obj);
            [U,S,V]=svd(eta);
            obj1=complexnumber2ellipticcomplex(obj,U);
            obj2=complexnumber2ellipticcomplex(obj,S);
            obj3=complexnumber2ellipticcomplex(obj,V);
        end

        %Returns pseudo inverse of the elliptic matrices.
        function obj=elliptic_n_pinv(obj)
            eta=ellipticcomplex2complexnumber(obj);
            cpinv=pinv(eta);
            obj=complexnumber2ellipticcomplex(obj,cpinv);
        end

        %Returns the rank of the elliptic matrices.
        function v=elliptic_n_rank(obj)
            eta=ellipticcomplex2complexnumber(obj);
            v=rank(eta);
        end

        %Returns the ith row elements of the elliptic matrices.
        function obj1=elliptic_n_i_th_erow(obj,i)
            eta=ellipticcomplex2complexnumber(obj);
            eta_A_ith_row=eta(i,:);
            obj1=complexnumber2ellipticcomplex(obj,eta_A_ith_row);
        end

        %Returns the ith column elements of the elliptic matrices.
        function obj1=elliptic_n_i_th_ecolumn(obj,i)
            eta=ellipticcomplex2complexnumber(obj);
            eta_A_ith_column=eta(:,i);
            obj1=complexnumber2ellipticcomplex(obj,eta_A_ith_column);
        end

        %Returns the ith row and jth column element of the elliptic
        %matrices.
        function obj1=elliptic_n_index(obj,i,j)
            eta=ellipticcomplex2complexnumber(obj);
            obj1=complexnumber2ellipticcomplex(obj,eta(i,j));
        end

        %Returns the least squares solution of the elliptic matrix equation
        %obj1*X=obj2 and the least square error.
        function [obj,error]=elliptic_n__lss(obj1,obj2)
            obj1_eta=ellipticcomplex2complexnumber(obj1);
            obj2_eta=ellipticcomplex2complexnumber(obj2);
            solution=pinv(obj1_eta)*obj2_eta;
            obj=complexnumber2ellipticcomplex(obj2,solution);
            error=e_norm(elliptic_n_subtraction(elliptic_n_product(obj1,obj),obj2));
        end

    end %methods

end %classdef