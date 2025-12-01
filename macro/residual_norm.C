#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "sdst.h" // 确保包含 DST 类的头文件

void normalize(int isotope, std::string input_path, std::string output_path, std::string fit_ref_file) {
    // 打开输入文件
    TFile *input = TFile::Open(input_path.c_str(), "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << input_path << std::endl;
        return;
    }

    // 获取输入树
    TTree *t = (TTree*)input->Get("tree");
    if (!t) {
        std::cerr << "Error: Cannot find tree 'tree' in input file" << std::endl;
        input->Close();
        return;
    }

    // 打开拟合参考文件
    TFile *fit_ref = TFile::Open(fit_ref_file.c_str(), "READ");
    if (!fit_ref || fit_ref->IsZombie()) {
        std::cerr << "Error: Cannot open fit reference file " << fit_ref_file << std::endl;
        input->Close();
        return;
    }

    // 预先加载所有拟合函数
    std::map<std::string, TF1*> mean_fits;
    std::map<std::string, TF1*> sigma_fits;
    
    for (int det = 0; det < 3; det++) {
        for (int layer = 1; layer < 8; layer++) {
            std::string mean_key = Form("f_mean_be%d_%d_%d", isotope, det, layer);
            std::string sigma_key = Form("f_sigma_be%d_%d_%d", isotope, det, layer);
            
            mean_fits[mean_key] = dynamic_cast<TF1*>(fit_ref->Get(mean_key.c_str()));
            sigma_fits[sigma_key] = dynamic_cast<TF1*>(fit_ref->Get(sigma_key.c_str()));
            
            if (!mean_fits[mean_key] || !sigma_fits[sigma_key]) {
                std::cerr << "Error: Missing fit function for " << mean_key << " or " << sigma_key << std::endl;
                input->Close();
                fit_ref->Close();
                return;
            }
            
            // 设置函数持久化，避免自动删除
            mean_fits[mean_key]->SetBit(kCanDelete, false);
            sigma_fits[sigma_key]->SetBit(kCanDelete, false);
        }
    }
    
    // 创建输出文件
    TFile *output = TFile::Open(output_path.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << output_path << std::endl;
        input->Close();
        fit_ref->Close();
        return;
    }
    
    // 创建克隆树（只复制结构，不复制数据）
    TTree *new_tree = t->CloneTree(0);
    
    // 设置DST对象
    DST *dst = new DST();
    dst->SetAddress(t);
    
    // 添加新分支
    Float_t y_residual[3][7]; // [detector][layer]
    new_tree->Branch("y_residual", y_residual, "y_residual[3][7]/F");
    
    // 处理所有条目
    const Long64_t nEntries = t->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        if (i % 10000 == 0) 
            std::cout << "Processing entry: " << i << "/" << nEntries << std::endl;
        
        t->GetEntry(i);
        
        for (int det = 0; det < 3; det++) {
            for (int layer = 1; layer < 8; layer++) {
                std::string mean_key = Form("f_mean_be%d_%d_%d", isotope, det, layer);
                std::string sigma_key = Form("f_sigma_be%d_%d_%d", isotope, det, layer);
                
                TF1 *mean_fit = mean_fits[mean_key];
                TF1 *sigma_fit = sigma_fits[sigma_key];
                
                double mean, sigma;
                if (det == 0) {
                    mean = mean_fit->Eval(dst->fEk.Ek[0]);
                    sigma = sigma_fit->Eval(dst->fEk.Ek[0]);
                } else {
                    mean = mean_fit->Eval(dst->fEk.Ek[1]);
                    sigma = sigma_fit->Eval(dst->fEk.Ek[1]);
                }
                
                if (sigma == 0) sigma = 1e-6;
                
                // 注意：layer 索引从1开始，但数组索引从0开始
                y_residual[det][layer-1] = ((dst->fTrack.cooy[layer] - dst->fTrack.yglob_JN[layer]) - mean) / sigma;
            }
        }
        
        new_tree->Fill();
    }
    
    // 写入并关闭文件
    output->cd();
    new_tree->Write();
    output->Close();
    input->Close();
    fit_ref->Close();
    
    // 清理内存
    delete dst;
    
    // 注意：不要删除TF1对象，它们属于fit_ref文件
    // 只需清除map（但不删除对象）
    mean_fits.clear();
    sigma_fits.clear();
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <isotope> <input_path> <output_path> <fit_ref_file>" << std::endl;
        return 1;
    }
    int isotope = atoi(argv[1]);
    std::string input_path = argv[2];
    std::string output_path = argv[3];
    std::string fit_ref_file = argv[4];

    normalize(isotope, input_path, output_path, fit_ref_file);
    std::cout << "Normalization completed successfully." << std::endl;
    return 0;
}