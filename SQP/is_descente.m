function [H, rho, check_reset, check_descente]=is_descente(grad_merite_x, H, rho, check_reset)
% La fonction teste si d est bien une direction de descente 
% et réinitialise H puis augmente rho si la réinitialisation
% de H n'a pas suffit.

        check_descente = (grad_merite_x < 0 );
        if ~check_descente && ~check_reset
            %Réinitialisation de la hessienne
            H=eye(size(H,1));
            check_reset=true;
        elseif ~(check_descente) && check_reset
            % Augmentation du paramètre de pénalisation
            rho=rho*10;
            assert(rho<10^(15))
        else
            % new_d est une direction de descente
            check_descente=true;
        end
end