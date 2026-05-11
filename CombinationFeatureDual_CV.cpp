#include "pch.h"
#include "CombinationFeatureDual_CV.h"

#include "DefineInspInfo.h"
#include "Recipe.h"
#include "inspectionType.h"
#include "ClassifierManagerDual.h"
#include "Common.h"

#if __has_include(<opencv2/core.hpp>)
#include <opencv2/core.hpp>
#define COMBINATION_CV_OPENCV_ENABLED 1
#else
#define COMBINATION_CV_OPENCV_ENABLED 0
#endif

namespace
{
	constexpr double kDefaultZero = 0.0;
}

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
	UNREFERENCED_PARAMETER(a_ctx);
	UNREFERENCED_PARAMETER(a_pattern_index);
	const Json::Value& layer_pattern_array = a_recipe[RCP_ROOT::VP][a_ctx.cam_index][RCP_VP::PTN];

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

#if COMBINATION_CV_OPENCV_ENABLED
			cv::Point2d defect_pt(col, row);
			double pseudo_feature = cv::norm(defect_pt) * 0.0;
#else
			double pseudo_feature = 0.0;
#endif
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitScore"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmitScore"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["BlackDOmitScore"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Omit_Ratio_64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmit_Ratio_64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Black_DOmit_Ratio_64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Pre64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMin_Pre64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMax_Pre64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GraySTDEV_Pre64"] = kDefaultZero;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["CV_OpenCV_Ready"] = 1;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["CV_OpenCV_Feature"] = pseudo_feature;
		}
	}
}
