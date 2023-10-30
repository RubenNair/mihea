# Initialize variables to store the sum and count of the columns you're interested in
column3_sum = 0
column5_sum = 0
count = 0

# Open the file and read its lines
with open("output/meeting_25_10/normalizedCVarsBasedOnNumSamples/statistics.txt", "r") as file:
    lines = file.readlines()

# Iterate through lines, skipping the first line
for line in lines[1:]:
    # Split the line into columns using space as a delimiter
    columns = line.split()
    if len(columns) >= 5:
        # Extract the third and fifth columns
        column3 = float(columns[2])
        column5 = float(columns[4].split("(")[1].split(",")[0])

        # Add the values to the sum
        column3_sum += column3
        column5_sum += column5

        # Increase the count
        count += 1

# Calculate the averages
average_column3 = column3_sum / count
average_column5 = column5_sum / count

print(f"Average of the third column: {average_column3}")
print(f"Average of the parsed fifth column: {average_column5}")
