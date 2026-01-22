void initial_setting_wri(
    const std::filesystem::path& dir, const std::filesystem::path& subdir, const std::ostringstream& file_name, 
    const unsigned int Lx, const unsigned int Ly, 
    std::vector<Par_i>& par_vec_i, std::vector<Par_f>& par_vec_f, std::vector<unsigned int>& t_vec, 
    const unsigned int* R1A, const unsigned int* R2A, const unsigned int* R3A, const unsigned int* R4A, const unsigned int* TA, 
    const unsigned int* R1B, const unsigned int* R2B, const unsigned int* R3B, const unsigned int* R4B, const unsigned int* TB 
) 
{
    H5::H5File file((dir / subdir / file_name.str()).string(), H5F_ACC_TRUNC); 
    // 
    H5::StrType str_type(H5::PredType::C_S1, 64); 
    str_type.setCset(H5T_CSET_UTF8); 
    str_type.setStrpad(H5T_STR_NULLTERM); 
    H5::CompType comp_type_i(sizeof(Par_i)); 
    H5::CompType comp_type_f(sizeof(Par_f)); 
    comp_type_i.insertMember("name", HOFFSET(Par_i, name), str_type); 
    comp_type_i.insertMember("value", HOFFSET(Par_i, value), H5::PredType::NATIVE_UINT); 
    comp_type_f.insertMember("name", HOFFSET(Par_f, name), str_type); 
    comp_type_f.insertMember("value", HOFFSET(Par_f, value), H5::PredType::NATIVE_DOUBLE); 
    hsize_t dims_i[1] = {par_vec_i.size()}; 
    H5::DataSpace dataspace_i(1, dims_i); 
    H5::DataSet dataset_i = file.createDataSet("parameters_int", comp_type_i, dataspace_i); 
    dataset_i.write(par_vec_i.data(), comp_type_i); 
    hsize_t dims_f[1] = {par_vec_f.size()}; 
    H5::DataSpace dataspace_f(1, dims_f); 
    H5::DataSet dataset_f = file.createDataSet("parameters_float", comp_type_f, dataspace_f); 
    dataset_f.write(par_vec_f.data(), comp_type_f); 
    hsize_t dims_time[1] = {t_vec.size()}; 
    H5::DataSpace dataspace_time(1, dims_time); 
    H5::DataSet dataset_time = file.createDataSet("time_vector", H5::PredType::NATIVE_UINT, dataspace_time); 
    dataset_time.write(t_vec.data(), H5::PredType::NATIVE_UINT); 
    hsize_t dims_res[3] = {1, Ly, Lx}; 
    // hsize_t dims_res[3] = {1, Lx, Ly}; 
    hsize_t maxdims_res[3] = {H5S_UNLIMITED, Ly, Lx}; 
    // hsize_t maxdims_res[3] = {H5S_UNLIMITED, Lx, Ly}; 
    hsize_t chunk_dims[3] = {1, Ly, Lx}; 
    // hsize_t chunk_dims[3] = {1, Lx, Ly}; 
    H5::DSetCreatPropList prop; 
    prop.setChunk(3, chunk_dims); 
    H5::DataSpace dataspace_res(3, dims_res, maxdims_res); 
    H5::DataSet res_R1A = file.createDataSet("results_R1A", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R2A = file.createDataSet("results_R2A", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R3A = file.createDataSet("results_R3A", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R4A = file.createDataSet("results_R4A", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_TA = file.createDataSet("results_TA", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R1B = file.createDataSet("results_R1B", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R2B = file.createDataSet("results_R2B", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R3B = file.createDataSet("results_R3B", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_R4B = file.createDataSet("results_R4B", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    H5::DataSet res_TB = file.createDataSet("results_TB", H5::PredType::NATIVE_UINT, dataspace_res, prop); 
    res_R1A.write(R1A, H5::PredType::NATIVE_UINT); 
    res_R2A.write(R2A, H5::PredType::NATIVE_UINT); 
    res_R3A.write(R3A, H5::PredType::NATIVE_UINT); 
    res_R4A.write(R4A, H5::PredType::NATIVE_UINT); 
    res_TA.write(TA, H5::PredType::NATIVE_UINT); 
    res_R1B.write(R1B, H5::PredType::NATIVE_UINT); 
    res_R2B.write(R2B, H5::PredType::NATIVE_UINT); 
    res_R3B.write(R3B, H5::PredType::NATIVE_UINT); 
    res_R4B.write(R4B, H5::PredType::NATIVE_UINT); 
    res_TB.write(TB, H5::PredType::NATIVE_UINT); 
    dataset_i.close(); 
    dataset_f.close(); 
    dataset_time.close(); 
    res_R1A.close(); 
    res_R2A.close(); 
    res_R3A.close(); 
    res_R4A.close(); 
    res_TA.close(); 
    res_R1B.close(); 
    res_R2B.close(); 
    res_R3B.close(); 
    res_R4B.close(); 
    res_TB.close(); 
    file.flush(H5F_SCOPE_GLOBAL); 
    file.close(); 
} 
void results_wri(
    const std::filesystem::path& dir, const std::filesystem::path& subdir, const std::ostringstream& file_name, 
    const unsigned int j, const unsigned int Lx, const unsigned int Ly,  
    const unsigned int* R1A, const unsigned int* R2A, const unsigned int* R3A, const unsigned int* R4A, const unsigned int* TA, 
    const unsigned int* R1B, const unsigned int* R2B, const unsigned int* R3B, const unsigned int* R4B, const unsigned int* TB 
) 
{
    H5::H5File file((dir / subdir / file_name.str()).string(), H5F_ACC_RDWR); 
    H5::DataSet res_R1A = file.openDataSet("results_R1A"); 
    H5::DataSet res_R2A = file.openDataSet("results_R2A"); 
    H5::DataSet res_R3A = file.openDataSet("results_R3A"); 
    H5::DataSet res_R4A = file.openDataSet("results_R4A"); 
    H5::DataSet res_TA = file.openDataSet("results_TA"); 
    H5::DataSet res_R1B = file.openDataSet("results_R1B"); 
    H5::DataSet res_R2B = file.openDataSet("results_R2B"); 
    H5::DataSet res_R3B = file.openDataSet("results_R3B"); 
    H5::DataSet res_R4B = file.openDataSet("results_R4B"); 
    H5::DataSet res_TB = file.openDataSet("results_TB"); 
    hsize_t new_size[3] = {j+1, Ly, Lx}; 
    // hsize_t new_size[3] = {j+1, Lx, Ly}; 
    res_R1A.extend(new_size); 
    res_R2A.extend(new_size); 
    res_R3A.extend(new_size); 
    res_R4A.extend(new_size); 
    res_TA.extend(new_size); 
    res_R1B.extend(new_size); 
    res_R2B.extend(new_size); 
    res_R3B.extend(new_size); 
    res_R4B.extend(new_size); 
    res_TB.extend(new_size); 
    H5::DataSpace filespace1A = res_R1A.getSpace(); 
    H5::DataSpace filespace2A = res_R2A.getSpace(); 
    H5::DataSpace filespace3A = res_R3A.getSpace(); 
    H5::DataSpace filespace4A = res_R4A.getSpace(); 
    H5::DataSpace filespaceTA = res_TA.getSpace(); 
    H5::DataSpace filespace1B = res_R1B.getSpace(); 
    H5::DataSpace filespace2B = res_R2B.getSpace(); 
    H5::DataSpace filespace3B = res_R3B.getSpace(); 
    H5::DataSpace filespace4B = res_R4B.getSpace(); 
    H5::DataSpace filespaceTB = res_TB.getSpace(); 
    hsize_t offset[3] = {j, 0, 0}; 
    hsize_t count[3] = {1, Ly, Lx}; 
    // hsize_t count[3] = {1, Lx, Ly}; 
    filespace1A.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace2A.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace3A.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace4A.selectHyperslab(H5S_SELECT_SET, count, offset);  
    filespaceTA.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace1B.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace2B.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace3B.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespace4B.selectHyperslab(H5S_SELECT_SET, count, offset); 
    filespaceTB.selectHyperslab(H5S_SELECT_SET, count, offset); 
    H5::DataSpace memspace(3, count); 
    res_R1A.write(R1A, H5::PredType::NATIVE_UINT, memspace, filespace1A); 
    res_R2A.write(R2A, H5::PredType::NATIVE_UINT, memspace, filespace2A); 
    res_R3A.write(R3A, H5::PredType::NATIVE_UINT, memspace, filespace3A); 
    res_R4A.write(R4A, H5::PredType::NATIVE_UINT, memspace, filespace4A); 
    res_TA.write(TA, H5::PredType::NATIVE_UINT, memspace, filespaceTA); 
    res_R1B.write(R1B, H5::PredType::NATIVE_UINT, memspace, filespace1B); 
    res_R2B.write(R2B, H5::PredType::NATIVE_UINT, memspace, filespace2B); 
    res_R3B.write(R3B, H5::PredType::NATIVE_UINT, memspace, filespace3B); 
    res_R4B.write(R4B, H5::PredType::NATIVE_UINT, memspace, filespace4B); 
    res_TB.write(TB, H5::PredType::NATIVE_UINT, memspace, filespaceTB); 
    memspace.close(); 
    filespace1A.close(); 
    filespace2A.close(); 
    filespace3A.close(); 
    filespace4A.close(); 
    filespaceTA.close(); 
    filespace1B.close(); 
    filespace2B.close(); 
    filespace3B.close(); 
    filespace4B.close(); 
    filespaceTB.close(); 
    res_R1A.close(); 
    res_R2A.close(); 
    res_R3A.close(); 
    res_R4A.close(); 
    res_TA.close(); 
    res_R1B.close(); 
    res_R2B.close(); 
    res_R3B.close(); 
    res_R4B.close(); 
    res_TB.close(); 
    file.flush(H5F_SCOPE_GLOBAL); 
    file.close(); 
}