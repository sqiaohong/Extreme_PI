function H_extract = H_matrix(his_time_len, his_fu_time_len)

    H_extract = zeros(his_time_len, his_fu_time_len);
    H_extract(:, 1:his_time_len) = eye(his_time_len);
end

