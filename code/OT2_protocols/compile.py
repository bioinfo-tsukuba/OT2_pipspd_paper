import sys


def output_common_imports():
    try:
        with open('commons.py', 'r') as file:
            for line in file:
                # end='' is used to avoid adding extra newlines
                print(line, end='')
    except FileNotFoundError:
        print("Error: File 'common_imports.py' not found.")
    except Exception as e:
        print(f"Error: An unexpected error occurred - {e}")


def output_after_end_of_imports(file_name):
    try:
        with open(file_name, 'r') as file:
            after_end_of_imports = False
            for line in file:
                if "# END_OF_IMPORTS" in line:
                    after_end_of_imports = True
                    continue
                if after_end_of_imports:
                    # end='' is used to avoid adding extra newlines
                    print(line, end='')
    except FileNotFoundError:
        print(f"Error: File '{file_name}' not found.")
    except Exception as e:
        print(f"Error: An unexpected error occurred - {e}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py file_name")
    else:
        file_name = sys.argv[1]
        output_common_imports()
        output_after_end_of_imports(file_name)
