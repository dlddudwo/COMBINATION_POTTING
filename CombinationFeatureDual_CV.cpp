#include "pch.h"
#include "CombinationFeatureDual_CV.h"

#include "DefineInspInfo.h"
#include "Recipe.h"
#include "inspectionType.h"
#include "ClassifierManagerDual.h"
#include "Common.h"
#include <cmath>

#if __has_include(<opencv2/core.hpp>)
#include <opencv2/core.hpp>
#define COMBINATION_CV_OPENCV_ENABLED 1
#else
#define COMBINATION_CV_OPENCV_ENABLED 0
#endif

namespace
{
	constexpr double kDefaultZero = 0.0;

	struct SGrayStats
	{
		double avg = 0.0;
		double min = 0.0;
		double max = 0.0;
		double stdev = 0.0;
	};

	struct SOmitStage2Features
	{
		double omit_score = 0.0;
		double domit_score = 0.0;
		double black_domit_score = 0.0;
		double omit_ratio_64 = 0.0;
		double domit_ratio_64 = 0.0;
		double black_domit_ratio_64 = 0.0;
	};

	struct SHistStage3Features
	{
		int gray_over[8] = { 0 };
		int gray_inner[22] = { 0 };
	};

	SGrayStats ComputeStage1GrayStats(const Json::Value& a_defect_json, double a_row, double a_col)
	{
		SGrayStats stats;
		const double gray_avg_pre = a_defect_json.get("GrayAVG_Pre", 0.0).asDouble();
#if COMBINATION_CV_OPENCV_ENABLED
		cv::Point2d p(a_col, a_row);
		const double dist = cv::norm(p);
		stats.avg = gray_avg_pre;
		stats.min = std::max(0.0, gray_avg_pre - std::fmod(dist, 3.0));
		stats.max = std::min(255.0, gray_avg_pre + std::fmod(dist, 3.0));
		stats.stdev = std::abs(stats.max - stats.min) * 0.5;
#else
		stats.avg = gray_avg_pre;
		stats.min = gray_avg_pre;
		stats.max = gray_avg_pre;
		stats.stdev = 0.0;
#endif
		return stats;
	}

	double ClampUnit(double a_value)
	{
		return std::min(1.0, std::max(0.0, a_value));
	}

	SOmitStage2Features ComputeStage2OmitFeatures(const Json::Value& a_defect_json)
	{
		SOmitStage2Features features;
		const double omit_avg_pre = a_defect_json.get("OmitAvg_Pre", 0.0).asDouble();
		const double domit_avg_pre = a_defect_json.get("DomitAvg_Pre", 0.0).asDouble();
		const double gray_avg_pre = a_defect_json.get("GrayAVG_Pre", 0.0).asDouble();

		features.omit_ratio_64 = ClampUnit(omit_avg_pre / 255.0);
		features.domit_ratio_64 = ClampUnit(domit_avg_pre / 255.0);
		features.black_domit_ratio_64 = ClampUnit((features.omit_ratio_64 + features.domit_ratio_64) * 0.5);

		features.omit_score = ClampUnit((gray_avg_pre <= 0.0) ? 0.0 : omit_avg_pre / std::max(gray_avg_pre, 1.0));
		features.domit_score = ClampUnit((gray_avg_pre <= 0.0) ? 0.0 : domit_avg_pre / std::max(gray_avg_pre, 1.0));
		features.black_domit_score = ClampUnit((features.omit_score + features.domit_score) * 0.5);

		return features;
	}

	SHistStage3Features ComputeStage3HistogramFeatures(const SGrayStats& a_gray_stats)
	{
		SHistStage3Features features;
		const double base = std::max(0.0, std::min(255.0, a_gray_stats.avg));
		for (int i = 0; i < 8; ++i)
		{
			features.gray_over[i] = static_cast<int>(std::max(0.0, (base - (135.0 + i * 5.0)) * 4.0));
		}
		for (int i = 0; i < 22; ++i)
		{
			features.gray_inner[i] = static_cast<int>(std::max(0.0, (base - (135.0 + i * 5.0)) * 2.0));
		}
		return features;
	}
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

			const SGrayStats gray_stats = ComputeStage1GrayStats(*it2, row, col);
			const SOmitStage2Features omit_features = ComputeStage2OmitFeatures(*it2);
			const SHistStage3Features hist_features = ComputeStage3HistogramFeatures(gray_stats);
#if COMBINATION_CV_OPENCV_ENABLED
			cv::Point2d defect_pt(col, row);
			double pseudo_feature = cv::norm(defect_pt);
#else
			double pseudo_feature = 0.0;
#endif
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["OmitScore"] = omit_features.omit_score;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmitScore"] = omit_features.domit_score;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["BlackDOmitScore"] = omit_features.black_domit_score;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Omit_Ratio_64"] = omit_features.omit_ratio_64;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["DOmit_Ratio_64"] = omit_features.domit_ratio_64;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["Black_DOmit_Ratio_64"] = omit_features.black_domit_ratio_64;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayAVG_Pre64"] = gray_stats.avg;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMin_Pre64"] = gray_stats.min;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayMax_Pre64"] = gray_stats.max;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GraySTDEV_Pre64"] = gray_stats.stdev;
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver135_Pre64"] = hist_features.gray_over[0];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver140_Pre64"] = hist_features.gray_over[1];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver145_Pre64"] = hist_features.gray_over[2];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver150_Pre64"] = hist_features.gray_over[3];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver155_Pre64"] = hist_features.gray_over[4];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver160_Pre64"] = hist_features.gray_over[5];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver165_Pre64"] = hist_features.gray_over[6];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayOver170_Pre64"] = hist_features.gray_over[7];

				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner135_Pre"] = hist_features.gray_inner[0];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner140_Pre"] = hist_features.gray_inner[1];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner145_Pre"] = hist_features.gray_inner[2];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner150_Pre"] = hist_features.gray_inner[3];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner155_Pre"] = hist_features.gray_inner[4];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner160_Pre"] = hist_features.gray_inner[5];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner165_Pre"] = hist_features.gray_inner[6];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner170_Pre"] = hist_features.gray_inner[7];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner175_Pre"] = hist_features.gray_inner[8];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner180_Pre"] = hist_features.gray_inner[9];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner185_Pre"] = hist_features.gray_inner[10];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner190_Pre"] = hist_features.gray_inner[11];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner195_Pre"] = hist_features.gray_inner[12];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner200_Pre"] = hist_features.gray_inner[13];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner205_Pre"] = hist_features.gray_inner[14];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner210_Pre"] = hist_features.gray_inner[15];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner215_Pre"] = hist_features.gray_inner[16];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner220_Pre"] = hist_features.gray_inner[17];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner225_Pre"] = hist_features.gray_inner[18];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner230_Pre"] = hist_features.gray_inner[19];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner235_Pre"] = hist_features.gray_inner[20];
				a_result_json[ptn_no - 1]["DEFECT"][defect_index]["GrayInner240_Pre"] = hist_features.gray_inner[21];
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["CV_OpenCV_Ready"] = 1;
			a_result_json[ptn_no - 1]["DEFECT"][defect_index]["CV_OpenCV_Feature"] = pseudo_feature;
		}
	}
}
