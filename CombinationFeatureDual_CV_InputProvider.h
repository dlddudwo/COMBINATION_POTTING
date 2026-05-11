#pragma once

#include <string>

#if __has_include(<opencv2/core.hpp>)
#include <opencv2/core.hpp>
#define COMBINATION_CV_HAS_OPENCV 1
#else
#define COMBINATION_CV_HAS_OPENCV 0
#endif

struct SCombinationDualCvInput
{
#if COMBINATION_CV_HAS_OPENCV
	cv::Mat img_pre;
	cv::Mat img_ori;
	cv::Mat mask_omit;
	cv::Mat mask_domit;
	cv::Mat mask_black_domit;
#endif
	bool valid = false;
};

class CCombinationFeatureDualCvInputProvider
{
public:
	SCombinationDualCvInput GetInput(const std::string& a_panel_id, int a_cam_index, int a_pattern_index) const;
};
