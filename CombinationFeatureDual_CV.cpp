#include "pch.h"
#include "CombinationFeatureDual_CV.h"

#include "DefineInspInfo.h"
#include "Recipe.h"
#include "inspectionType.h"
#include "SharedMemoryCollection.h"
#include "FeatureName.h"
#include "Common.h"
#include "HalconMath.h"
#include "ClassifierManagerDual.h"

Json::Value CCombinationFeatureDual_CV::Run(Json::Value a_recipe, Json::Value& a_result_json)
{
	CTimeChecker tmc(__FUNCTION__);

	const SRunContext ctx = BuildRunContext(a_recipe, a_result_json);
	const SPatternTypeIndex pattern_index = FindPatternTypeIndex(a_recipe, ctx.cam_index);

	ProcessPatterns(a_recipe, a_result_json, ctx, pattern_index);

	return a_result_json;
}

CCombinationFeatureDual_CV::SRunContext CCombinationFeatureDual_CV::BuildRunContext(const Json::Value& a_recipe, const Json::Value& a_result_json) const
{
	SRunContext ctx;

	ctx.cam_num = a_result_json[0]["VpNo"].asInt();
	ctx.cam_index = max(ctx.cam_num - 1, 0);

	const Json::Value& layer_insp_info = a_recipe[RCP_ROOT::INSP_INFO];
	ctx.pattern_count = layer_insp_info[RCP_INSP_INFO::PTN_COUNT].asInt();
	ctx.panel_id = layer_insp_info[RCP_INSP_INFO::PANEL_NAME].asString();

	return ctx;
}

CCombinationFeatureDual_CV::SPatternTypeIndex CCombinationFeatureDual_CV::FindPatternTypeIndex(const Json::Value& a_recipe, int a_cam_index) const
{
	SPatternTypeIndex pattern_index;

	const Json::Value& layer_pattern_array = a_recipe[RCP_ROOT::VP][a_cam_index][RCP_VP::PTN];
	for (int i = 0; i < layer_pattern_array.size(); ++i)
	{
		const PATTERN_TYPE pattern_type = static_cast<PATTERN_TYPE>(layer_pattern_array[i][RCP_PTN::PTN_TYPE].asInt());
		if (pattern_type == PATTERN_TYPE::OMIT) pattern_index.omit = i;
		if (pattern_type == PATTERN_TYPE::DOMIT) pattern_index.domit = i;
		if (pattern_type == PATTERN_TYPE::BLACK_DOMIT) pattern_index.black_domit = i;
	}

	return pattern_index;
}

void CCombinationFeatureDual_CV::ProcessPatterns(const Json::Value& a_recipe, Json::Value& a_result_json, const SRunContext& a_ctx, const SPatternTypeIndex& a_pattern_index)
{
	const Json::Value& layer_pattern_array = a_recipe[RCP_ROOT::VP][a_ctx.cam_index][RCP_VP::PTN];

	CHalconMath halcon;
	auto& shared_memory = CSharedMemoryCollection::GetInstance();
	int error = 0;
	const int fixed_ref_inspection_type = 0;
	HObject region_inspection = shared_memory.GetInspectionRegionHalcon(a_ctx.panel_id, a_ctx.cam_index, INSPECTION_REGION_TYPE::NORMAL, error);

	std::unique_ptr<HObject[]> arr_ho_original(new HObject[a_ctx.pattern_count]);
	std::unique_ptr<HObject[]> arr_ho_pre_processing(new HObject[a_ctx.pattern_count]);
	for (size_t i = 0; i < a_ctx.pattern_count; ++i)
	{
		arr_ho_original[i] = shared_memory.GetImageHalcon(a_ctx.panel_id, a_ctx.cam_index, IMAGE_TYPE::ORIGINAL, i, 0);
		arr_ho_pre_processing[i] = shared_memory.GetImageHalcon(a_ctx.panel_id, a_ctx.cam_index, IMAGE_TYPE::PRE_PROCESSING, i, fixed_ref_inspection_type);
	}

	HObject ho_binarized_omit = shared_memory.GetImageHalcon(a_ctx.panel_id, a_ctx.cam_index, IMAGE_TYPE::PRE_PROCESSING, a_pattern_index.omit, 0);
	HObject ho_binarized_domit = shared_memory.GetImageHalcon(a_ctx.panel_id, a_ctx.cam_index, IMAGE_TYPE::PRE_PROCESSING, a_pattern_index.domit, 0);
	HObject ho_binarized_black_domit = shared_memory.GetImageHalcon(a_ctx.panel_id, a_ctx.cam_index, IMAGE_TYPE::PRE_PROCESSING, a_pattern_index.black_domit, 0);
	UNREFERENCED_PARAMETER(ho_binarized_omit);
	UNREFERENCED_PARAMETER(ho_binarized_domit);
	UNREFERENCED_PARAMETER(ho_binarized_black_domit);

	HObject ho_omit = shared_memory.GetMask(a_ctx.panel_id, a_ctx.cam_index, MASK_TYPE::OMIT, error);
	HObject ho_domit = shared_memory.GetMask(a_ctx.panel_id, a_ctx.cam_index, MASK_TYPE::DOMIT, error);
	HObject ho_black_domit = shared_memory.GetMask(a_ctx.panel_id, a_ctx.cam_index, MASK_TYPE::BLACK_DOMIT, error);
	Union1(ho_omit, &ho_omit);
	Union1(ho_domit, &ho_domit);
	Union1(ho_black_domit, &ho_black_domit);
	UNREFERENCED_PARAMETER(halcon);
	UNREFERENCED_PARAMETER(region_inspection);
	UNREFERENCED_PARAMETER(ho_omit);
	UNREFERENCED_PARAMETER(ho_domit);
	UNREFERENCED_PARAMETER(ho_black_domit);

	for (Json::ValueIterator it = a_result_json.begin(); it != a_result_json.end(); ++it)
	{
		const int pattern_index = (*it)["PTNNo"].asInt() - 1;
		const PATTERN_TYPE pattern_type = static_cast<PATTERN_TYPE>(layer_pattern_array[pattern_index][RCP_PTN::PTN_TYPE].asInt());
		if (pattern_type == PATTERN_TYPE::OMIT || pattern_type == PATTERN_TYPE::DOMIT || pattern_type == PATTERN_TYPE::BLACK_DOMIT) continue;

		const int ptn_no = (*it)["PTNNo"].asInt();
		Json::Value ptn_json = (*it)["DEFECT"];
		int defect_index = -1;

		for (Json::ValueIterator it2 = ptn_json.begin(); it2 != ptn_json.end(); ++it2)
		{
			if ((*it2).isMember(CDualClassifierFeature::Classify_Group) == false) continue;
			defect_index++;

			double row = (*it2)["Row"].asDouble();
			double col = (*it2)["Column"].asDouble();
			const int resize_ratio = (*it2)["Resize_Ratio"].asDouble();
			if (resize_ratio > 0)
			{
				row *= resize_ratio;
				col *= resize_ratio;
			}

			const int blob_region_index = (*it2)["Region_Index"].asInt();
			HObject extract_region;
			HObject extract_region_origin;
			HObject region64;
			if (blob_region_index < 0)
			{
				GenRegionPoints(&extract_region, row, col);
				DilationRectangle1(extract_region, &extract_region, 3, 3);
			}
			else
			{
				extract_region = shared_memory.GetBlobRegion(a_ctx.panel_id, a_ctx.cam_index, ptn_no - 1, blob_region_index);
			}

			GenRegionPoints(&region64, row, col);
			DilationRectangle1(region64, &region64, 64, 64);
			Intersection(region64, region_inspection, &region64);

			if (resize_ratio > 0) extract_region_origin = halcon.ResizeRegion(extract_region, resize_ratio);
			else extract_region_origin = extract_region;

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
			if (halcon.ValidHRegion(ho_domit) == true)
			{
				Intersection(ho_domit, extract_region_origin, &ho_intersection);
				AreaCenter(ho_intersection, &hv_area, &hv_row_omit, &hv_col_omit);
				double score = hv_area.D() / hv_area_region.D();
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmitScore"] = score;
			}
			if (halcon.ValidHRegion(ho_black_domit) == true)
			{
				Intersection(ho_black_domit, extract_region_origin, &ho_intersection);
				AreaCenter(ho_intersection, &hv_area, &hv_row_omit, &hv_col_omit);
				double score = hv_area.D() / hv_area_region.D();
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["BlackDOmitScore"] = score;
			}

			HTuple hv_pre_average64_omit, hv_pre_dsd64_omit;
			HTuple hv_pre_average64_domit, hv_pre_dsd64_domit;
			HTuple hv_pre_average64_black_domit, hv_pre_dsd64_black_domit;
			Intensity(region64, ho_binarized_omit, &hv_pre_average64_omit, &hv_pre_dsd64_omit);
			double omit_ratio = hv_pre_average64_omit.D() / 255;
			Intensity(region64, ho_binarized_domit, &hv_pre_average64_domit, &hv_pre_dsd64_domit);
			double domit_ratio = hv_pre_average64_domit.D() / 255;
			Intensity(region64, ho_binarized_black_domit, &hv_pre_average64_black_domit, &hv_pre_dsd64_black_domit);
			double black_domit_ratio = hv_pre_average64_black_domit.D() / 255;

			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Omit_Ratio_64"] = omit_ratio;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmit_Ratio_64"] = domit_ratio;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Black_DOmit_Ratio_64"] = black_domit_ratio;

			HTuple hv_pre_omit_average, hv_pre_Domit_average = 0.0;
			HTuple hv_ori_omit_average, hv_ori_Domit_average = 0.0;
			if (halcon.ValidHRegion(ho_omit) == true)
			{
				Intensity(extract_region, arr_ho_original[a_pattern_index.omit], &hv_ori_omit_average, NULL);
				Intensity(extract_region, arr_ho_pre_processing[a_pattern_index.omit], &hv_pre_omit_average, NULL);
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitAvg_Pre"] = (double)hv_pre_omit_average;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitAvg_Ori"] = (double)hv_ori_omit_average;
			}
			if (a_pattern_index.domit != 0)
			{
				Intensity(extract_region, arr_ho_original[a_pattern_index.domit], &hv_ori_Domit_average, NULL);
				Intensity(extract_region, arr_ho_pre_processing[a_pattern_index.domit], &hv_pre_Domit_average, NULL);
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DomitAvg_Pre"] = (double)hv_pre_Domit_average;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DomitAvg_Ori"] = (double)hv_ori_Domit_average;
			}

			HTuple hv_PreAverage64, hv_PredSd64, hv_PreMin64, hv_PreMax64, hv_PreRange64;
			MinMaxGray(region64, arr_ho_pre_processing[ptn_no - 1], 0, &hv_PreMin64, &hv_PreMax64, &hv_PreRange64);
			Intensity(region64, arr_ho_pre_processing[ptn_no - 1], &hv_PreAverage64, &hv_PredSd64);
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Pre64"] = (double)hv_PreAverage64;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMin_Pre64"] = (double)hv_PreMin64;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMax_Pre64"] = (double)hv_PreMax64;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GraySTDEV_Pre64"] = (double)hv_PredSd64;
		}
	}
}
