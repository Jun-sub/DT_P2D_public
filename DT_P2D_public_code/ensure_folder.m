function ensure_folder(folderPath)
    if ~isfolder(folderPath)
        mkdir(folderPath);
    end
end
