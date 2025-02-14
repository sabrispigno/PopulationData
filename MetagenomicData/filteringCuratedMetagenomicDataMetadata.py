import pandas as pd

class MetadataFilter:
    def __init__(self, metadata):
        """
        Initialize the MetadataFilter with a DataFrame.

        :param metadata: pd.DataFrame containing the metadata to filter
        """
        self.metadata = metadata

    def filter_data(self, study_name=None, subject_id=None, body_sites=None, antibiotics_current_use=None, 
                    diets=None,
                    study_condition=None, diseases=None, age=None, infant_age=None, age_categories=None, 
                    genders=None, countries=None, non_westernized=None, sequencing_platform=None, 
                    DNA_extraction_kit=None, PMID=None, number_reads=None, number_bases=None, 
                    minimum_read_length=None, median_read_length=None, NCBI_accession=None):
        """
        Filter the metadata based on the specified criteria.
        You should provide the criteria as a list of values to include.
        """
        filtered_data = self.metadata

        # Apply filters based on the presence of arguments
        if study_name:
            filtered_data = filtered_data[filtered_data["study_name"].isin(study_name)]
        if subject_id:
            filtered_data = filtered_data[filtered_data["subject_id"].isin(subject_id)]
        if body_sites:
            filtered_data = filtered_data[filtered_data["body_site"].isin(body_sites)]
        if antibiotics_current_use:
            filtered_data = filtered_data[filtered_data["antibiotics_current_use"] != antibiotics_current_use]
        if diets:
            filtered_data = filtered_data[filtered_data["diet"].isin(diets)]
        if study_condition:
            filtered_data = filtered_data[filtered_data["study_condition"].isin(study_condition)]
        if diseases:
            filtered_data = filtered_data[filtered_data["disease"].isin(diseases)]
        if age:
            filtered_data = filtered_data[filtered_data["age"].isin(age)]
        if infant_age:
            filtered_data = filtered_data[filtered_data["infant_age"].isin(infant_age)]
        if age_categories:
            filtered_data = filtered_data[filtered_data["age_category"].isin(age_categories)]
        if genders:
            filtered_data = filtered_data[filtered_data["gender"].isin(genders)]
        if countries:
            filtered_data = filtered_data[filtered_data["country"].isin(countries)]
        if non_westernized:
            filtered_data = filtered_data[filtered_data["non_westernized"].isin(non_westernized)]
        if sequencing_platform:
            filtered_data = filtered_data[filtered_data["sequencing_platform"].isin(sequencing_platform)]
        if DNA_extraction_kit:
            filtered_data = filtered_data[filtered_data["DNA_extraction_kit"].isin(DNA_extraction_kit)]
        if PMID:
            filtered_data = filtered_data[filtered_data["PMID"].isin(PMID)]
        if number_reads:
            filtered_data = filtered_data[filtered_data["number_reads"].isin(number_reads)]
        if number_bases:
            filtered_data = filtered_data[filtered_data["number_bases"].isin(number_bases)]
        if minimum_read_length:
            filtered_data = filtered_data[filtered_data["minimum_read_length"].isin(minimum_read_length)]
        if median_read_length:
            filtered_data = filtered_data[filtered_data["median_read_length"].isin(median_read_length)]
        if NCBI_accession:
            filtered_data = filtered_data[filtered_data["NCBI_accession"].isin(NCBI_accession)]

        return filtered_data