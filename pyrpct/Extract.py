# import packages
import os
import copy
from Feature import feature_kpssm, feature_dtpssm, extract_reduce
from Visual import easy_time, visual_aa
from Read import read_raac, read_pssm

extract_path = os.path.dirname(__file__)


# 获取raac
def extract_raa(raa, outfolder, now_path):
    # 获取氨基酸约化密码表
    raa_path = os.path.join(extract_path, 'raacDB')
    raa_file = os.path.join(raa_path, raa)
    if raa in os.listdir(raa_path):
        raacode = read_raac(raa_file)
        if outfolder not in os.listdir(now_path):
            outfolder = os.path.join(now_path, outfolder)
            os.makedirs(outfolder)
        else:
            outfolder = os.path.join(now_path, outfolder)
    else:
        with open(raa_file, 'w') as f:
            f.write('type 1 size ' + str(len(raa.split('-'))) + ' ' + raa)
        raacode = read_raac(raa_file)
        if outfolder not in os.listdir(now_path):
            outfolder = now_path
        else:
            outfolder = now_path
        os.remove(raa_file)
    return raacode, outfolder


# 提取矩阵特征
def extract_features(pssm_matrixes, pssm_aaid):
    all_features = []
    start_e = 0
    for i in range(len(pssm_matrixes)):
        start_e += 1
        easy_time(start_e, len(pssm_matrixes))
        each_matrix = pssm_matrixes[i]
        matrix_400 = []
        aa_index = visual_aa()
        for aa in aa_index:
            aa_score = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0]
            for j in range(len(each_matrix)):
                line = each_matrix[j]
                if pssm_aaid[i][j] == aa:
                    for k in range(len(line)):
                        aa_score[k] = aa_score[k] + line[k]
            matrix_400.append(aa_score)
        all_features.append(matrix_400)
    return all_features


# 矩阵转置
def extract_transform(data):
    new_data = []
    for i in data:
        new_data += i
    return new_data


# 归一化
def extract_scale(data):
    new_data = []
    min_x = min(data)
    max_x = max(data)
    for i in data:
        if (max_x - min_x) != 0:
            new_data.append(round((i - min_x) / (max_x - min_x), 4))
        else:
            new_data.append(0)
    return new_data


# min max
def extract_minmax(data):
    new_data = copy.deepcopy(data)
    for i in range(len(data)):
        for j in range(len(data[i])):
            min_x = min(data[i][j])
            max_x = max(data[i][j])
            for k in range(len(data[i][j])):
                if (max_x - min_x) != 0:
                    new_data[i][j][k] = round((data[i][j][k] - min_x) / (max_x - min_x), 4)
                else:
                    new_data[i][j][k] = 0
    return new_data


# 保存特征
def extract_save(raa_features, kpssm_features, dtpssm_features, outfolder, raa_list):
    start_e = 0
    for k in range(len(raa_features)):
        start_e += 1
        easy_time(start_e, len(raa_features))
        eachraa = raa_features[k]
        out_file = ''
        for i in range(len(eachraa)):
            eachfile = eachraa[i]
            # eachaac = aac_features[i]
            eachkpssm = kpssm_features[i][k]
            # eachkmer = psekraac_features[i][k]
            # eachsaac = saac_features[i]
            # eachsw = sw_features[i][k]
            eachdtpssm = dtpssm_features[i][k]
            type_m = eachfile[0]
            data_m = eachfile[1]
            data_m = extract_transform(data_m) + eachkpssm + eachdtpssm
            # 归一化
            data_m = extract_scale(data_m)
            mid_file = type_m
            for j in range(len(data_m)):
                mid_file += ' ' + str(j + 1) + ':' + str(data_m[j])
            out_file += mid_file + '\n'
        path = os.path.join(outfolder, raa_list[k] + '_rpct.fs')
        with open(path, 'w') as f2:
            f2.write(out_file)
            f2.close()


# extract main
def extract_main(positive, negative, outfolder, raa, lmda, now_path):
    # 获取raac
    raacode, outfolder = extract_raa(raa, outfolder, now_path)
    # 处理地址
    positive = os.path.join(os.path.join(now_path, 'PSSMs'), positive)
    negative = os.path.join(os.path.join(now_path, 'PSSMs'), negative)
    # positive
    pssm_matrixes, pssm_aaid, pssm_type = read_pssm(positive, [], [], [], '0')
    # negative
    pssm_matrixes, pssm_aaid, pssm_type = read_pssm(negative, pssm_matrixes, pssm_aaid, pssm_type, '1')
    # PSSM特征提取
    pssm_features = extract_features(pssm_matrixes, pssm_aaid)
    pssm_features = extract_minmax(pssm_features)
    pssm_matrixes = extract_minmax(pssm_matrixes)
    # Sequence PseKRAAC
    # k, g, m = 2, 1, 1
    # psekraac_features = feature_psekraac(raacode, pssm_aaid, k, g, m)
    # AAC特征提取
    # aac_features = feature_aac(pssm_aaid)
    # SAAC特征提取
    # saac_features = feature_saac(pssm_aaid)
    # kpssm特征提取
    kpssm_features = feature_kpssm(pssm_matrixes, raacode)
    # psepssm特征提取
    dtpssm_features = feature_dtpssm(pssm_matrixes, raacode, 3)
    # PSSM滑窗
    # sw_features = feature_sw(pssm_matrixes, pssm_aaid, raacode, lmda)
    # 矩阵约化
    raa_features = extract_reduce(pssm_features, raacode, pssm_type)
    # 生成特征文件
    extract_save(raa_features, kpssm_features, dtpssm_features, outfolder, raacode[1])
