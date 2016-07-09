#include <gtest/gtest.h>
#include "./test_models.cpp"
#include "./test_loadvcf.cpp"
#include "./test_liabilitymodel.cpp"

int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
