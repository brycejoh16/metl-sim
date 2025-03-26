import pandas as pd
import unittest
import os
from merge_rosetta_energies import join_dataframes


class TestJoinDataframes(unittest.TestCase):
    def setUp(self):
        """Set up test environment."""
        self.test_output_dir = "test_output"
        os.makedirs(self.test_output_dir, exist_ok=True)
        self.github_commit = "abc123"

        # Create sample DataFrames and save as CSV
        self.df1 = pd.DataFrame({
            "variant": ["v1", "v2", "v3"],
            "energy": [1.0, 2.0, 3.0]
        })
        self.df2 = pd.DataFrame({
            "variant": ["v1", "v3", "v4"],
            "energy": [1.5, 2.5, 3.5]
        })
        self.df3 = pd.DataFrame({
            "variant": ["v1", "v2", "v3"],
            "constant_col": [10, 10, 10]  # Zero variance column
        })

        self.paths = [os.path.join(self.test_output_dir, f"df_{i}.csv") for i in range(1, 4)]
        self.suffixes = ["a", "b", "c"]

        # Save CSVs
        self.df1.to_csv(self.paths[0], index=False)
        self.df2.to_csv(self.paths[1], index=False)
        self.df3.to_csv(self.paths[2], index=False)

    def tearDown(self):
        """Clean up after test."""
        for path in self.paths:
            os.remove(path)
        if os.path.exists(self.test_output_dir):
            for file in os.listdir(self.test_output_dir):
                os.remove(os.path.join(self.test_output_dir, file))
            os.rmdir(self.test_output_dir)

    def test_join_dataframes(self):
        """Test the join_dataframes function."""
        log_path = os.path.join(self.test_output_dir, "merge_log.txt")
        output_path = os.path.join(self.test_output_dir, "energies_df.csv")

        # Run the function
        join_dataframes(self.paths, self.suffixes, self.test_output_dir, self.github_commit)

        # Check if the output file is created
        self.assertTrue(os.path.exists(output_path))

        # Load the merged DataFrame
        merged_df = pd.read_csv(output_path)

        # Check if zero variance columns are removed
        self.assertNotIn("constant_col_c", merged_df.columns)

        # Check if the correct variants are merged
        expected_variants = ["v1", "v3"]
        self.assertCountEqual(merged_df["variant"].tolist(), expected_variants)

        # Check if the log file is created and has content
        self.assertTrue(os.path.exists(log_path))
        with open(log_path, "r") as log_file:
            log_content = log_file.read()
            self.assertIn("GitHub Commit: abc123", log_content)
            self.assertIn("Loaded", log_content)
            self.assertIn("Removing columns with zero variance", log_content)
            self.assertIn("Final merged dataframe saved", log_content)

if __name__ == '__main__':
    unittest.main()