% 读取Excel文件
file_path = 'C:\Users\Admin\Desktop\Data of dravel curve_Origin_Data.xlsx';
data = readtable(file_path);

% 获取数据大小
[num_rows, num_cols] = size(data);

% 设置每组的行数
group_size = 2;

% 计算需要多少组
num_groups = ceil(num_rows / group_size);

% 创建新表格存储平均值
avg_table = table();

% 对每一列进行处理
for col = 1:num_cols
    % 提取当前列数据（假设是数值数据）
    column_data = table2array(data(:, col));
    
    % 用于存储每组的平均值
    group_means = zeros(num_groups, 1);
    
    % 对每组数据计算平均值
    for group = 1:num_groups
        % 计算当前组的起始和结束索引
        start_idx = (group - 1) * group_size + 1;
        end_idx = min(group * group_size, num_rows);
        
        % 提取当前组的数据
        group_data = column_data(start_idx:end_idx);
        
        % 计算平均值
        group_means(group) = mean(group_data);
    end
    
    % 将平均值添加到新表格
    col_name = data.Properties.VariableNames{col};
    avg_table.(col_name) = group_means;
end

% 保存为新的Excel文件
output_file = 'C:\Users\Admin\Desktop\Data of dravel curve_Origin_Data_2.xlsx';
writetable(avg_table, output_file);

% 显示结果
disp('处理完成！');
disp(['原始数据行数: ', num2str(num_rows)]);
disp(['平均值数据行数: ', num2str(num_groups)]);
disp(['输出文件: ', output_file]);

% 显示前几行结果预览
disp('前5行平均值结果:');
disp(head(avg_table, 5));