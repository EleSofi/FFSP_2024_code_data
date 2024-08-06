function main()
    % Automatically determine the path to the current script's directory
    current_script_dir = fileparts(mfilename('fullpath'));

    % Determine the project root directory
    try
        project_root = find_project_root(current_script_dir, 'FFSP_code_revised_2024');
        disp(['Project root found: ' project_root]);
    catch ME
        disp(ME.message);
        project_root = [];
    end

    if ~isempty(project_root)
        output_file = fullfile(current_script_dir, 'project_structure.md');
        generate_directory_structure(project_root, output_file);
    else
        disp('Project root not found. Project structure not generated.');
    end
end

function generate_directory_structure(root_dir, output_file)
    fid = fopen(output_file, 'w');
    if fid == -1
        error('Cannot open output file.');
    end
    
    fprintf(fid, '# Project Structure\n\n');
    write_structure(fid, root_dir, root_dir, 0);
    
    fclose(fid);
    disp(['Project structure saved to ' output_file]);
end

function write_structure(fid, root_dir, current_dir, level)
    indent = repmat(' ', 1, 4 * level);
    sub_indent = repmat(' ', 1, 4 * (level + 1));
    
    fprintf(fid, '%s- %s/\n', indent, get_relative_path(root_dir, current_dir));
    
    files = dir(current_dir);
    for i = 1:length(files)
        if files(i).isdir
            if ~strcmp(files(i).name, '.') && ~strcmp(files(i).name, '..')
                write_structure(fid, root_dir, fullfile(current_dir, files(i).name), level + 1);
            end
        else
            fprintf(fid, '%s- %s\n', sub_indent, files(i).name);
        end
    end
end

function relative_path = get_relative_path(root_dir, current_dir)
    relative_path = strrep(current_dir, [root_dir, filesep], '');
end

function project_root = find_project_root(starting_dir, root_marker)
    current_dir = starting_dir;
    while true
        if any(strcmp({dir(current_dir).name}, root_marker))
            project_root = fullfile(current_dir, root_marker);
            return;
        end
        parent_dir = fileparts(current_dir);
        if strcmp(parent_dir, current_dir)
            error('Project root directory ''%s'' not found.', root_marker);
        end
        current_dir = parent_dir;
    end
end

% Call the main function
main();
