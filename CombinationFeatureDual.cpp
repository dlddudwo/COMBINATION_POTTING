#include "pch.h"
#include "CombinationFeatureDual.h"
#include "SharedMemoryCollection.h"
#include "FeatureName.h"
#include "DefineInspInfo.h"
#include "Recipe.h"
#include "Common.h"
#include "HalconMath.h"
#include "inspectionType.h"
#include "ClassifierManagerDual.h"

Json::Value CCombinationFeatureDual::Run(Json::Value a_recipe, Json::Value& a_result_json)
{
	CTimeChecker tmc(__FUNCTION__);

	const int cam_num = a_result_json[0]["VpNo"].asInt();
	const int cam_index = max(cam_num - 1, 0);

	Json::Value& layer_insp_info = a_recipe[RCP_ROOT::INSP_INFO];
	auto& layer_cam = a_recipe[RCP_ROOT::VP][cam_index];
	auto& layer_pattern_array = layer_cam[RCP_VP::PTN];

	std::string model_name = layer_insp_info[RCP_INSP_INFO::MODEL_NAME].asString();
	std::string panel_id = layer_insp_info[RCP_INSP_INFO::PANEL_NAME].asString();

	CHalconMath halcon;
	auto& config = CConfig::GetInstance();
	auto& shared_memory = CSharedMemoryCollection::GetInstance();
	try
	{
		int vpNo = cam_num;
		//string panel_id = panel_name;
		int pattern_count = layer_insp_info[RCP_INSP_INFO::PTN_COUNT].asInt();

		vector<bool> vec_priority_check;

		string judge = "";

		int error = 0;
		const int fixed_ref_inspection_type = 0;
		HObject region_inspection = shared_memory.GetInspectionRegionHalcon(panel_id, cam_index, INSPECTION_REGION_TYPE::NORMAL, error);
		bool use_max_contrast_70 = false;

		std::unique_ptr<HObject[]> arr_ho_original(new HObject[pattern_count]);
		std::unique_ptr<HObject[]> arr_ho_pre_processing(new HObject[pattern_count]);
		std::unique_ptr<HObject[]> arr_ho_fft(new HObject[pattern_count]);
		for (size_t i = 0; i < pattern_count; ++i)
		{
			arr_ho_original[i] = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::ORIGINAL, i, 0);
			arr_ho_pre_processing[i] = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::PRE_PROCESSING, i, fixed_ref_inspection_type);
			if (use_max_contrast_70) arr_ho_fft[i] = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::FFT_DENOISE, i, 0);
		}

		// Omit 유형별 index 탐색
		int omit_index = 0, domit_index = 0, black_domit_index = 0;
		for (int i = 0; i < layer_pattern_array.size(); i++)
		{
			PATTERN_TYPE pattern_type = static_cast<PATTERN_TYPE>(layer_pattern_array[i][RCP_PTN::PTN_TYPE].asInt());
			if (pattern_type == PATTERN_TYPE::OMIT)			omit_index = i;
			if (pattern_type == PATTERN_TYPE::DOMIT)		domit_index = i;
			if (pattern_type == PATTERN_TYPE::BLACK_DOMIT)	black_domit_index = i;
		}

		HObject ho_binarized_omit = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::PRE_PROCESSING, omit_index, 0);
		HObject ho_binarized_domit = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::PRE_PROCESSING, domit_index, 0);
		HObject ho_binarized_black_domit = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::PRE_PROCESSING, black_domit_index, 0);

		HObject ho_omit = shared_memory.GetMask(panel_id, cam_index, MASK_TYPE::OMIT, error);
		HObject ho_domit = shared_memory.GetMask(panel_id, cam_index, MASK_TYPE::DOMIT, error);
		HObject ho_black_domit = shared_memory.GetMask(panel_id, cam_index, MASK_TYPE::BLACK_DOMIT, error);

		Union1(ho_omit, &ho_omit);
		Union1(ho_domit, &ho_domit);
		Union1(ho_black_domit, &ho_black_domit);
		// HObject HResizeRgn;
		// GenEmptyObj(&HResizeRgn);

		for (Json::ValueIterator it = a_result_json.begin(); it != a_result_json.end(); it++)// PTN
		{
			int pattern_index = (*it)["PTNNo"].asInt() - 1;
			PATTERN_TYPE pattern_type = static_cast<PATTERN_TYPE>(layer_pattern_array[pattern_index][RCP_PTN::PTN_TYPE].asInt());
			if (pattern_type == PATTERN_TYPE::OMIT || pattern_type == PATTERN_TYPE::DOMIT || pattern_type == PATTERN_TYPE::BLACK_DOMIT) continue; // Omit PTN들 제외
			int ptn_no = (*it)["PTNNo"].asInt();

			Json::Value ptn_json = (*it)["DEFECT"];
			int defect_index = -1;

			for (Json::ValueIterator it2 = ptn_json.begin(); it2 != ptn_json.end(); it2++)
			{
				if ((*it2).isMember(CDualClassifierFeature::Classify_Group) == false) continue;

				defect_index++;

				double connection = (*it2)[FEATURE_NAME::Connection].asInt();

				double avg_gray = (*it2)["GrayAVG_Pre"].asDouble();

				double row = (*it2)["Row"].asDouble();
				double col = (*it2)["Column"].asDouble();

				int resize_ratio = (*it2)["Resize_Ratio"].asDouble();
				if (resize_ratio > 0)
				{
					row *= resize_ratio;
					col *= resize_ratio;
				}

				int blob_region_index = (*it2)["Region_Index"].asInt();
				int ptn_type = (*it2)["PTNType"].asInt();
				int inspection_type_index = (*it2)["inspection_type_index"].asInt();

				HObject extract_region;
				HObject extract_region_origin;
				HObject region64;
				HObject region64_bg;
				HObject extract_region64;
				int inspection_region_index = ptn_no - 1;

				GenRegionPoints(&region64, row, col);
				DilationRectangle1(region64, &region64, 64, 64); // 일단 64로 하드코딩
				Intersection(region64, region_inspection, &region64); // 외곽 처리를 위해 Intersection 함

				if (blob_region_index < 0)
				{
					GenRegionPoints(&extract_region, row, col);
					DilationRectangle1(extract_region, &extract_region, 3, 3);
				}
				else
				{
					extract_region = shared_memory.GetBlobRegion(panel_id, cam_index, ptn_no - 1, blob_region_index);
				}
				//-----------------------------------------------------------------------------------------
				if (resize_ratio > 0)
					extract_region_origin = halcon.ResizeRegion(extract_region, resize_ratio);
				else
					extract_region_origin = extract_region;

				// 0. Omit Score 계산
				// 1. Omit 
				// 2. Omit 밝기 연산 
				HObject ho_intersection;
				HTuple hv_row_omit, hv_col_omit, hv_area, hv_row_region, hv_col_region, hv_area_region;
				AreaCenter(extract_region_origin, &hv_area_region, &hv_row_region, &hv_col_region);

				if (halcon.ValidHRegion(ho_omit) == true)
				{
					Intersection(ho_omit, extract_region_origin, &ho_intersection);
					AreaCenter(ho_intersection, &hv_area, &hv_row_omit, &hv_col_omit);
					double score = hv_area.D() / hv_area_region.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitScore"] = score;
				}
				// 2. DOmit
				if (halcon.ValidHRegion(ho_domit) == true)
				{
					Intersection(ho_domit, extract_region_origin, &ho_intersection);
					AreaCenter(ho_intersection, &hv_area, &hv_row_omit, &hv_col_omit);
					double score = hv_area.D() / hv_area_region.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmitScore"] = score;
				}
				// 3. Black+DOmit
				if (halcon.ValidHRegion(ho_black_domit) == true)
				{
					Intersection(ho_black_domit, extract_region_origin, &ho_intersection);
					AreaCenter(ho_intersection, &hv_area, &hv_row_omit, &hv_col_omit);
					double score = hv_area.D() / hv_area_region.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["BlackDOmitScore"] = score;
				}
				//-----------------------------------------------------------------------------------------
				// 0. 64x64 영역 Omit 추출
				HTuple hv_pre_average64_omit, hv_pre_dsd64_omit;
				HTuple hv_pre_average64_domit, hv_pre_dsd64_domit;
				HTuple hv_pre_average64_black_domit, hv_pre_dsd64_black_domit;

				//HObject ho_binarized_omit = shared_memory.GetImageHalcon(panel_name, cam_index, IMAGE_TYPE::PRE_PROCESSING, 0, 0);
				Intensity(region64, ho_binarized_omit, &hv_pre_average64_omit, &hv_pre_dsd64_omit);
				double omit_ratio = hv_pre_average64_omit.D() / 255;
				Intensity(region64, ho_binarized_domit, &hv_pre_average64_domit, &hv_pre_dsd64_domit);
				double domit_ratio = hv_pre_average64_domit.D() / 255;
				Intensity(region64, ho_binarized_black_domit, &hv_pre_average64_black_domit, &hv_pre_dsd64_black_domit);
				double black_domit_ratio = hv_pre_average64_black_domit.D() / 255;

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Omit_Ratio_64"] = omit_ratio;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmit_Ratio_64"] = domit_ratio;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Black_DOmit_Ratio_64"] = black_domit_ratio;
				//-----------------------------------------------------------------------------------------
				//0. Omit 영역내 밝기 값 추출   불필요 중복코드 확인여부 필요(25.10.21)
				HTuple hv_pre_omit_average, hv_pre_Domit_average = 0.0;
				HTuple hv_ori_omit_average, hv_ori_Domit_average = 0.0;

				if (halcon.ValidHRegion(ho_omit) == true)
				{
					//Intersection(ho_omit, extract_region, &ho_intersection);
					Intensity(extract_region, arr_ho_original[omit_index], &hv_ori_omit_average, NULL);
					Intensity(extract_region, arr_ho_pre_processing[omit_index], &hv_pre_omit_average, NULL);
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitAvg_Pre"] = (double)hv_pre_omit_average;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitAvg_Ori"] = (double)hv_ori_omit_average;
				}
				// 2. DOmit
				if (domit_index != 0)
				{
					//Intersection(ho_domit, extract_region, &ho_intersection);
					Intensity(extract_region, arr_ho_original[domit_index], &hv_ori_Domit_average, NULL);
					Intensity(extract_region, arr_ho_pre_processing[domit_index], &hv_pre_Domit_average, NULL);
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DomitAvg_Pre"] = (double)hv_pre_Domit_average;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DomitAvg_Ori"] = (double)hv_ori_Domit_average;
				}

				// -----------------------------------------------------------------------------------------
				// 0. 64x64 영역 Feature 추출
				HTuple hv_PreAverage64, hv_PredSd64, hv_PreMin64, hv_PreMax64, hv_PreRange64;
				HTuple hv_PreAverage64_bg, hv_PredSd64_bg;

				MinMaxGray(region64, arr_ho_pre_processing[ptn_no - 1], 0, &hv_PreMin64, &hv_PreMax64, &hv_PreRange64);
				Intensity(region64, arr_ho_pre_processing[ptn_no - 1], &hv_PreAverage64, &hv_PredSd64);

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Pre64"] = (double)hv_PreAverage64;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMin_Pre64"] = (double)hv_PreMin64;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMax_Pre64"] = (double)hv_PreMax64;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GraySTDEV_Pre64"] = (double)hv_PredSd64;

				// 8은 ptn개수 아님 135 ~ 170까지 8 step이라 그냥 8로 하드코딩해버림


				HTuple hv_AbsHisto_GrayOver, hv_ReHist_GrayOver, hv_AbsHisto_GrayInner, hv_ReHisto_GrayInner;
				HTuple hv_select;
				//GrayOver
				int arrGrayOver[8];
				int hist_GrayOver[256] = { 0 };
				int cum_step_GrayOver[256] = { 0 };
				//GrayInner
				int hist_GrayInner[256] = { 0 };
				int cum_step_GrayInner[256] = { 0 };
				int arrGrayInner[22];

				int startGray = 135;
				int step = 5;

				GrayHisto(region64, arr_ho_pre_processing[ptn_no - 1], &hv_AbsHisto_GrayOver, &hv_ReHist_GrayOver);
				GrayHisto(extract_region, arr_ho_pre_processing[ptn_no - 1], &hv_AbsHisto_GrayInner, &hv_ReHisto_GrayInner);

				for (int i = 0; i < 256; ++i)
				{
					hist_GrayOver[i] = hv_AbsHisto_GrayOver[i].I();
					hist_GrayInner[i] = hv_AbsHisto_GrayInner[i].I();
					// step별 누적합 생성
					if (i == 0)
					{
						cum_step_GrayOver[i] = hist_GrayOver[i];
						cum_step_GrayInner[i] = hist_GrayInner[i];
					}
					else
					{
						cum_step_GrayOver[i] = cum_step_GrayOver[i - 1] + hist_GrayOver[i];
						cum_step_GrayInner[i] = cum_step_GrayInner[i - 1] + hist_GrayInner[i];
					}
				}

				// 구간 합 계산
				for (int i = 0; i < 8; ++i)
				{
					arrGrayOver[i] = cum_step_GrayOver[255] - cum_step_GrayOver[startGray - 1];
					startGray += step;
				}
				startGray = 135;
				// 구간 합 계산
				for (int i = 0; i < 22; ++i)
				{
					arrGrayInner[i] = cum_step_GrayInner[255] - cum_step_GrayInner[startGray - 1];
					startGray += step;
				}
				//GrayOver
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver135_Pre64"] = arrGrayOver[0];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver140_Pre64"] = arrGrayOver[1];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver145_Pre64"] = arrGrayOver[2];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver150_Pre64"] = arrGrayOver[3];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver155_Pre64"] = arrGrayOver[4];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver160_Pre64"] = arrGrayOver[5];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver165_Pre64"] = arrGrayOver[6];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver170_Pre64"] = arrGrayOver[7];

				//GrayInner
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner135_Pre"] = arrGrayInner[0];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner140_Pre"] = arrGrayInner[1];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner145_Pre"] = arrGrayInner[2];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner150_Pre"] = arrGrayInner[3];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner155_Pre"] = arrGrayInner[4];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner160_Pre"] = arrGrayInner[5];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner165_Pre"] = arrGrayInner[6];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner170_Pre"] = arrGrayInner[7];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner175_Pre"] = arrGrayInner[8];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner180_Pre"] = arrGrayInner[9];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner185_Pre"] = arrGrayInner[10];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner190_Pre"] = arrGrayInner[11];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner195_Pre"] = arrGrayInner[12];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner200_Pre"] = arrGrayInner[13];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner205_Pre"] = arrGrayInner[14];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner210_Pre"] = arrGrayInner[15];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner215_Pre"] = arrGrayInner[16];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner220_Pre"] = arrGrayInner[17];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner225_Pre"] = arrGrayInner[18];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner230_Pre"] = arrGrayInner[19];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner235_Pre"] = arrGrayInner[20];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner240_Pre"] = arrGrayInner[21];

				// Org 비율 Feature
				HObject ho_Bin64Adaptive, ho_Bin64AdaptiveReduced;
				HTuple hv_Area_AdaptH, hv_Row_AdaptH, hv_Col_AdaptH;

				Intersection(region64, extract_region, &extract_region64);
				Difference(region64, extract_region64, &region64_bg);
				Intensity(region64_bg, arr_ho_original[ptn_no - 1], &hv_PreAverage64_bg, &hv_PredSd64_bg);
				double oriAVGGray = a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Org"].asDouble();
				double bgAVGGray = hv_PreAverage64_bg.D();
				double grayAVG_Rate_Org = oriAVGGray / bgAVGGray;

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Rate_Org"] = grayAVG_Rate_Org;

				double offTh = 6;
				int thLow = (int)(bgAVGGray + offTh) > 254 ? 254 : (int)(bgAVGGray + offTh);

				ReduceDomain(arr_ho_original[ptn_no - 1], region64, &ho_Bin64AdaptiveReduced);
				Threshold(ho_Bin64AdaptiveReduced, &ho_Bin64Adaptive, thLow, 255);
				AreaCenter(ho_Bin64Adaptive, &hv_Area_AdaptH, &hv_Row_AdaptH, &hv_Col_AdaptH);

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Area_Adaptive_H"] = hv_Area_AdaptH.L();

				//-----------------------------------------------------------------------------------------
				// 1. PTN 조합 Feature 추출

				for (int i = 0; i < pattern_count; i++)
				{
					int ptn_num = i + 1;
					HTuple hv_Average, hv_dSd, hv_Min, hv_Max, hv_Range;
					HTuple hv_PreAverage, hv_PredSd, hv_PreMin, hv_PreMax, hv_PreRange;

					MinMaxGray(extract_region_origin, arr_ho_original[i], 0, &hv_Min, &hv_Max, &hv_Range);
					Intensity(extract_region_origin, arr_ho_original[i], &hv_Average, &hv_dSd);

					MinMaxGray(extract_region_origin, arr_ho_pre_processing[i], 0, &hv_PreMin, &hv_PreMax, &hv_PreRange);
					Intensity(extract_region_origin, arr_ho_pre_processing[i], &hv_PreAverage, &hv_PredSd);

					string strMinOri = "MinGrayOri_PTN" + std::to_string(ptn_num);
					string strMaxOri = "MaxGrayOri_PTN" + std::to_string(ptn_num);
					string strAvgOri = "AvgGrayOri_PTN" + std::to_string(ptn_num);
					string strSTDEVOri = "STDEVGrayOri_PTN" + std::to_string(ptn_num);

					string strMinPre = "MinGrayPre_PTN" + std::to_string(ptn_num);
					string strMaxPre = "MaxGrayPre_PTN" + std::to_string(ptn_num);
					string strAvgPre = "AvgGrayPre_PTN" + std::to_string(ptn_num);
					string strSTDEVPre = "STDEVGrayPre_PTN" + std::to_string(ptn_num);

					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strMinOri] = (double)hv_Min;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strMaxOri] = (double)hv_Max;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strAvgOri] = (double)hv_Average;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strSTDEVOri] = (double)hv_dSd;

					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strMinPre] = (double)hv_PreMin;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strMaxPre] = (double)hv_PreMax;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strAvgPre] = (double)hv_PreAverage;
					a_result_json[ptn_no - 1]["DEFECT"][defect_index][strSTDEVPre] = (double)hv_PredSd;
				}
				//휘점계 후보점 중심 기준 인접 RGB Pixel 찾기
				std::string pattern_name = layer_pattern_array[pattern_index][RCP_PTN::PTN_NAME].asString();
				const int pixel_direction = layer_insp_info[RCP_INSP_INFO::PIXEL_DIRECTION].asInt();
				HTuple hv_meanstemp, hv_meantemp, hv_Deviationtemp, hv_MaxVal, hv_MaxIndex;
				HObject ho_linetemp;
				int arr_row[3], hRowInt, hColInt;
				if (pattern_name == "BLACK" || pattern_name == "BLACK_FL")
				{
					HTuple HRowInt = hv_row_region.TupleRound();
					HTuple HColInt = hv_col_region.TupleRound();
					for (int i = 1; i < 4; i++)//R,G,B PTN Index
					{
						for (int j = 0; j < 5; j++)
						{
							if (pixel_direction == 1)
							{
								hRowInt = HRowInt.I();
								hColInt = HColInt.I();
								GenRegionLine(&ho_linetemp, hRowInt + j - 2, hColInt - 32, hRowInt + j - 2, hColInt + 32);
							}
							else
							{
								hRowInt = HColInt.I();
								hColInt = HRowInt.I();
								GenRegionLine(&ho_linetemp, hColInt - 32, hRowInt + j - 2, hColInt + 32, hRowInt + j - 2);
							}
							Intersection(ho_linetemp, region_inspection, &ho_linetemp);
							Intensity(ho_linetemp, arr_ho_original[i], &hv_meantemp, &hv_Deviationtemp);
							hv_meanstemp[j] = hv_meantemp;
						}
						TupleMax(hv_meanstemp, &hv_MaxVal);
						TupleFind(hv_meanstemp, hv_MaxVal, &hv_MaxIndex);

						int maxIdx = hv_MaxIndex[0].I();
						int max_row = hRowInt + maxIdx - 1;
						arr_row[i - 1] = max_row;
					}

					int closestIndex = 0;
					int minDiff = std::abs(arr_row[0] - hRowInt);

					for (int i = 1; i < 3; i++)
					{
						int currentDiff = std::abs(arr_row[i] - hRowInt);
						if (currentDiff < minDiff)
						{
							minDiff = currentDiff;
							closestIndex = i;
						}
					}
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["RGBINFO_PTNType"] = closestIndex + 1;
				}
				else
				{
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["RGBINFO_PTNType"] = 0;
				}
				//-----------------------------------------------------------------------------------------
				// 2. 2nd Threshold
				HObject imgReduced, imgTh2ndBlack, imgTh2ndWhite, imgConcatTh;
				HTuple hv_Area, hv_row, hv_col, hv_length;
				HTuple hv_Average_2nd, hv_dSd_2nd, hv_Min_2nd, hv_Max_2nd, hv_Range_2nd, hv_Value;

				int insp_type_index = a_result_json[ptn_no - 1]["DEFECT"][defect_index]["inspection_type_index"].asInt();
				Json::Value insp_type_recipe = a_recipe[RCP_ROOT::VP][cam_index][RCP_VP::PTN][ptn_type][RCP_PTN::INSP_TYPE][insp_type_index];

				HObject ho_preprocessing_insp_type = shared_memory.GetImageHalcon(panel_id, cam_index, IMAGE_TYPE::PRE_PROCESSING, ptn_no - 1, insp_type_index);

				ReduceDomain(ho_preprocessing_insp_type, extract_region, &imgReduced);
				Threshold(imgReduced, &imgTh2ndBlack, 0, insp_type_recipe[RCP_INSP_TYPE::BTH_2ND].asInt());
				Threshold(imgReduced, &imgTh2ndWhite, insp_type_recipe[RCP_INSP_TYPE::WTH_2ND].asInt(), 255);
				ConcatObj(imgTh2ndBlack, imgTh2ndWhite, &imgConcatTh);
				Union1(imgConcatTh, &imgConcatTh);
				AreaCenter(imgConcatTh, &hv_Area, &hv_row, &hv_col);

				MinMaxGray(imgConcatTh, ho_preprocessing_insp_type, 0, &hv_Min_2nd, &hv_Max_2nd, &hv_Range_2nd);
				Intensity(imgConcatTh, ho_preprocessing_insp_type, &hv_Average_2nd, &hv_dSd_2nd);

				RegionFeatures(imgConcatTh, (HTuple("width").Append("height")), &hv_Value);

				int area2nd = 0;
				TupleLength(hv_Area, &hv_length);
				if (hv_length.I() > 0)
				{
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Area_2nd"] = hv_Area.L();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Pre_2nd"] = hv_Average_2nd.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GraySTDEV_Pre_2nd"] = hv_dSd_2nd.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMax_Pre_2nd"] = hv_Max_2nd.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMin_Pre_2nd"] = hv_Min_2nd.D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Width_2nd"] = hv_Value[0].D();
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Height_2nd"] = hv_Value[1].D();
				}
//#pragma omp for
				for (int i = 1; i < 6; i++)
				{
					if (insp_type_recipe[RCP_INSP_TYPE::BTH].asInt() != 0)
					{
						Threshold(imgReduced, &imgTh2ndBlack, 0, insp_type_recipe[RCP_INSP_TYPE::BTH].asInt() - i);
						/*					Threshold(imgReduced, &imgTh2ndWhite, insp_type_recipe[RCP_INSP_TYPE::WTH].asInt() + i, 255);
											ConcatObj(imgTh2ndBlack, imgTh2ndWhite, &imgConcatTh);
											Union1(imgConcatTh, &imgConcatTh);*/
						AreaCenter(imgTh2ndBlack, &hv_Area, &hv_row, &hv_col);
						a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Area-" + std::to_string(i)] = hv_Area.L();
					}
					else
						a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Area-" + std::to_string(i)] = -1;


				}
				//-----------------------------------------------------------------------------------------
				// 3. 4 Direction Feature 추출
				int pitch_x = insp_type_recipe[RCP_INSP_TYPE::PITCH_X].asInt();
				int pitch_y = insp_type_recipe[RCP_INSP_TYPE::PITCH_Y].asInt();
				HObject movedRegion[4];
				MoveRegion(extract_region_origin, &movedRegion[0], 0, -pitch_x);
				MoveRegion(extract_region_origin, &movedRegion[1], 0, pitch_x);
				MoveRegion(extract_region_origin, &movedRegion[2], pitch_y, 0);
				MoveRegion(extract_region_origin, &movedRegion[3], -pitch_y, 0);

				double totalGrayVal = 0;
				int count = 0;
//#pragma omp for
				for (int i = 0; i < 4; i++)
				{
					HObject ho_intersectRegion;
					HTuple hv_equal, hv_mean, hv_deviation, hv_meanSrc, hv_deviationSrc;
					Intersection(arr_ho_original[ptn_no - 1], movedRegion[i], &ho_intersectRegion);
					TestEqualRegion(movedRegion[i], ho_intersectRegion, &hv_equal);
					if (hv_equal.L() == 0) continue;

					Intensity(movedRegion[i], arr_ho_original[ptn_no - 1], &hv_mean, &hv_deviation);
					Intensity(extract_region_origin, arr_ho_original[ptn_no - 1], &hv_meanSrc, &hv_deviationSrc);

					totalGrayVal += abs(hv_mean.D() - hv_meanSrc.D());
					count++;
				}
				if (count > 0) totalGrayVal /= count;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayLevel_NealPxl_Org"] = totalGrayVal;
				//-----------------------------------------------------------------------------------------
				// 4. 상/하위 50% 정보 추출
				HTuple hv_regionRows, hv_regionCols, hv_regionGrayVals, hv_regionGrayValsSorted, hv_regionGrayValsSelected;
				HTuple hv_meanOriLT, hv_meanOriGT;
				hv_length.Clear();
				// Ori
				GetRegionPoints(extract_region_origin, &hv_regionRows, &hv_regionCols);
				GetGrayval(arr_ho_original[ptn_no - 1], hv_regionRows, hv_regionCols, &hv_regionGrayVals);
				TupleSort(hv_regionGrayVals, &hv_regionGrayValsSorted);
				TupleLength(hv_regionGrayValsSorted, &hv_length);
				hv_regionGrayValsSelected = hv_regionGrayValsSorted.TupleSelectRange(hv_length.I() / 2, hv_length.I() - 1);
				hv_meanOriGT = hv_regionGrayValsSelected.TupleMean();
				hv_regionGrayValsSelected = hv_regionGrayValsSorted.TupleSelectRange(0, hv_length.I() / 2 - 1);
				hv_meanOriLT = hv_regionGrayValsSelected.TupleMean();
				// Pre
				HTuple hv_regionRows2, hv_regionCols2, hv_regionGrayVals2, hv_regionGrayValsSorted2, hv_regionGrayValsSelected2;
				HTuple hv_meanPreLT, hv_meanPreGT;
				hv_length.Clear();

				GetRegionPoints(extract_region_origin, &hv_regionRows2, &hv_regionCols2);
				GetGrayval(arr_ho_pre_processing[ptn_no - 1], hv_regionRows2, hv_regionCols2, &hv_regionGrayVals2);
				TupleSort(hv_regionGrayVals2, &hv_regionGrayValsSorted2);
				TupleLength(hv_regionGrayValsSorted2, &hv_length);
				hv_regionGrayValsSelected2 = hv_regionGrayValsSorted2.TupleSelectRange(hv_length.I() / 2, hv_length.I() - 1);
				hv_meanPreGT = hv_regionGrayValsSelected2.TupleMean();
				hv_regionGrayValsSelected2 = hv_regionGrayValsSorted2.TupleSelectRange(0, hv_length.I() / 2 - 1);
				hv_meanPreLT = hv_regionGrayValsSelected2.TupleMean();

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["AvgGrayLevel_Org_50_LT"] = hv_meanOriLT.D();
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["AvgGrayLevel_Org_50_GT"] = hv_meanOriGT.D();
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["AvgGrayLevel_Pre_50_LT"] = hv_meanPreLT.D();
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["AvgGrayLevel_Pre_50_GT"] = hv_meanPreGT.D();

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayLevel_Rate"] = hv_meanOriLT.D() / totalGrayVal * 100;
				//-----------------------------------------------------------------------------------------
				// Max Contrast 70 추출
				bool use_max_contrast_70 = false;
				if (use_max_contrast_70)
				{

					const double contrast_ratio = 0.7;
					HTuple hv_gray_arr, hv_gray_arr_sorted, hv_pre_gray_avg, hv_pre_gray_dsd, hv_tuple_length, hv_tuple_selected, hv_tuple_mean;
					Intensity(extract_region, arr_ho_fft[ptn_no - 1], &hv_pre_gray_avg, &hv_pre_gray_dsd);
					GetGrayval(arr_ho_original[ptn_no - 1], hv_regionRows, hv_regionCols, &hv_gray_arr);
					TupleSort(hv_gray_arr, &hv_gray_arr_sorted);
					TupleLength(hv_gray_arr_sorted, &hv_tuple_length);

					if (hv_pre_gray_avg.D() < 127) // 암점
						hv_tuple_selected = hv_gray_arr_sorted.TupleSelectRange((int)(hv_tuple_length.I() * (1 - contrast_ratio)), hv_tuple_length.I() - 1);
					else // 휘점
						hv_tuple_selected = hv_gray_arr_sorted.TupleSelectRange(0, (int)(hv_tuple_length.I() * contrast_ratio) - 1);

					hv_tuple_mean = hv_tuple_selected.TupleMean();

					double dilation_px = 5.5;
					HObject ho_local_background, ho_local_background_diff;
					DilationCircle(extract_region, &ho_local_background, dilation_px);
					Difference(ho_local_background, extract_region, &ho_local_background_diff);
					Intensity(ho_local_background_diff, arr_ho_fft[ptn_no - 1], &hv_pre_gray_avg, &hv_pre_gray_dsd);
					// LocalBackground - grayDefectMax70Mean
					a_result_json[ptn_no - 1]["DEFECT"][defect_index]["MaxContrast70"] = hv_pre_gray_avg.D() - hv_tuple_mean.D();
				}
			}
		}
	}
	catch (Json::Exception& e)
	{
		const char* err_msg = e.what();
		spdlog::error("Json Exception Caught:{}", err_msg);
	}

	return a_result_json;

}
