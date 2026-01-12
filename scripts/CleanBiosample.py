import xml.etree.ElementTree as ET
import json

def parse_xml_to_dict(xml_file):
    """
    Parse XML file and create a dictionary mapping Name to [Name, HarmonizedName, Synonym1, Synonym2, ...]
    
    Args:
        xml_file: XML file path
    
    Returns:
        dict: Dictionary mapping Name to a deduplicated list containing Name, HarmonizedName and Synonyms
    """
    # Parse XML file
    tree = ET.parse(xml_file)
    root = tree.getroot()
    
    # Create dictionary to store mappings
    name_synonym_dict = {}
    
    # Iterate through all Attribute elements
    for attribute in root.findall('Attribute'):
        # Get Name as dictionary key
        name_elem = attribute.find('Name')
        if name_elem is not None and name_elem.text:
            name = name_elem.text
            
            # Create set to store Name, HarmonizedName and all Synonyms (automatically deduplicates)
            synonyms_set = set()
            
            # Add Name itself
            synonyms_set.add(name)
            
            # Add HarmonizedName
            harmonized_name_elem = attribute.find('HarmonizedName')
            if harmonized_name_elem is not None and harmonized_name_elem.text:
                synonyms_set.add(harmonized_name_elem.text)
            
            # Add all Synonyms
            for synonym_elem in attribute.findall('Synonym'):
                if synonym_elem.text:
                    synonyms_set.add(synonym_elem.text)
            
            # Convert set to list and store in dictionary
            name_synonym_dict[name] = list(synonyms_set)
    
    return name_synonym_dict

def save_dict_to_json(data_dict, output_file):
    """
    Save dictionary to JSON file
    
    Args:
        data_dict: Dictionary to save
        output_file: Output JSON file path
    """
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(data_dict, f, ensure_ascii=False, indent=2)
    print(f"Dictionary successfully saved to {output_file}")
    print(f"Processed {len(data_dict)} attributes")

if __name__ == "__main__":
    # Input file path
    xml_file = '/Users/a1-6/Project-helsinki/Meta2Data/scripts/Data/download.xml'
    # Output file path
    output_file = '/Users/a1-6/Project-helsinki/Meta2Data/scripts/Data/NCBI_Biosample.json'
    
    # Parse XML and create dictionary
    print("Parsing XML file...")
    name_synonym_dict = parse_xml_to_dict(xml_file)
    
    # Save to JSON file
    print("Saving to JSON file...")
    save_dict_to_json(name_synonym_dict, output_file)
    
    # Display first 3 entries as examples
    print("\nFirst 3 entries as examples:")
    for i, (key, values) in enumerate(list(name_synonym_dict.items())[:3]):
        print(f"{i+1}. '{key}': {values}\n")