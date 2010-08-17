add_test("Warnings" cmake -P Warnings.cmake)
add_test("Violations" cmake -P Violation.cmake)
add_subdirectory(gun)
add_subdirectory(ppbench)