% Foo Bar
clear;
clc;

%This is a MATLAB program that will execute the c++ executable that
%runs the list viterbi algorithm. The program takes the specified field
%values and creates a command-line argument that executes
%ListDecodingProject.exe with the required field entries. The results are
%returned in a labeled text file.

%The following variables must be instantiated:
%	int codeRateDenominator;
%	std::vector<int> trellisPolynomial;
%	int constraintLength;
%	int crcPolynomial;
%	int crcDegree;
%	int numOfSNRsComputed;
%	std::vector<double> snrDBList;
%	int numOfTrials;
%	int lengthOfSequence;
%	int maxListSize;
%	int minNumberOfErrors;

%Structure of the text file:
% First Line:
%(1)Number of Trials 
%(2)Length of input Sequence 
%(3)SNR 
%(4)Number of Bit Errors 
%(5)Number of Frame Errors 
%(6) Number of sequences that checked the CRC 
%(7)Largest List Size Value

trellisPolynomial = [133, 171];
codeRateDenominator = length(trellisPolynomial);
constraintLength = length(flip(de2bi(oct2dec(trellisPolynomial(1)))));
crcPolynomial = 3051;
crcPolynomialBinary = flip(de2bi(oct2dec(crcPolynomial)));
crcDegree = length(crcPolynomialBinary) - 1;
snrDBList = [-4:0.2:4];
numOfSNRsComputed = length(snrDBList);
numOfTrials = 10^4;
k = 2^4;
maxListSize = 10^6;
minNumberOfErrors = 250;
inputSequenceLength = 0;
tailBitingFlag = 0;
hardDecodingFlag = 0;

punctureFlag = 0;
punctureRepetitionLength= 7;
punctureNodeLocations = [0,9];
punctureCount = length(punctureNodeLocations);
if ~(punctureFlag == 0)
    additionalPunctureLocations = [];
else 
    additionalPunctureLocations = [];
end
additionalPunctureCount = length(additionalPunctureLocations);

%-------------------------------------------------------
% a = 0:0;
% b = dec2bin(a);
% inputSequenceList = [];
% lengthOfSequence = 8;
% 
% for iter = 1:length(b)
%    temp = [];
%    difference = lengthOfSequence-length(b(iter,:));
%    for iter2 = 1:difference
%       temp(iter2) = 0; 
%    end
%    for iter2 = 1:length(b(iter,:))
%       temp(iter2+difference) = str2num(b(iter,iter2)); 
%    end
%    inputSequenceList = [inputSequenceList;temp];
% end
%-------------------------------------------------------

inputSequenceList = 0;

for iter100 = 1:length(inputSequenceList(:,1))
    inputSequence = inputSequenceList(iter100,:);

    if (~(inputSequenceLength == 0))
        k = inputSequenceLength;
    end

    %Creating command line argument for the 
    filePath = '.\x64\Release\ListDecodingProject.exe';

    for snrIter = 1:length(snrDBList)

    %Creating array for simulation data
    sim_data = [];
    
        for iter0 = 1:length(maxListSize)
            commandLineArgs = '';
            commandLineArgs = commandLineArgs + string(codeRateDenominator);

            for iter = 1:length(trellisPolynomial)
                commandLineArgs = commandLineArgs + " " + string(trellisPolynomial(iter));
            end

            commandLineArgs = commandLineArgs + " " + string(constraintLength);
            commandLineArgs = commandLineArgs + " " + string(crcPolynomial);
            commandLineArgs = commandLineArgs + " " + string(crcDegree);
            commandLineArgs = commandLineArgs + " " + string(1);
            commandLineArgs = commandLineArgs + " " + string(snrDBList(snrIter));
            commandLineArgs = commandLineArgs + " " + string(numOfTrials);
            commandLineArgs = commandLineArgs + " " + string(k);
            commandLineArgs = commandLineArgs + " " + string(maxListSize(iter0));
            commandLineArgs = commandLineArgs + " " + string(minNumberOfErrors);
            commandLineArgs = commandLineArgs + " " + string(tailBitingFlag);
            commandLineArgs = commandLineArgs + " " + string(hardDecodingFlag);
            commandLineArgs = commandLineArgs + " " + string(punctureFlag);

            if (punctureFlag == 1)
                commandLineArgs = commandLineArgs + " " + string(punctureCount);
                commandLineArgs = commandLineArgs + " " + string(punctureRepetitionLength);
                for iter2 = 1:length(punctureNodeLocations)
                    commandLineArgs = commandLineArgs + " " + string(punctureNodeLocations(iter2));
                end
            end
            
            commandLineArgs = commandLineArgs + " " + string(additionalPunctureCount);
            
           if ~(additionalPunctureCount == 0)
                for iter2 = 1:length(additionalPunctureLocations)
                    commandLineArgs = commandLineArgs + " " + string(additionalPunctureLocations(iter2));
                end
            end
            
            if (~(inputSequenceLength == 0))
                commandLineArgs = commandLineArgs + " " + string(inputSequenceLength);
                for iterSequence = 1:inputSequenceLength
                    commandLineArgs = commandLineArgs + " " + string(inputSequence(iterSequence));
                end
            end

            command = filePath + " " + commandLineArgs;

            tic;
            system(command);
            toc;

            %Reading the values from the textfiles produced by C++ program
            sim_data = struct('num_of_trials',0,'max_list_size',0,'k',0, 'n',0, 'snr',0.0,'num_bit_errors',0, ...
            'num_frame_errors',0,'num_crc_check',0,'largest_list_depth',0,'listDepth',[],...
            'num_ue', [], 'num_nack', [], 'num_fe', [], 'total_insertions', 0, 'time_trellis', 0, 'time_traceback', 0);
            textFileName = "T";
            for polyIter = 1:length(trellisPolynomial)
                textFileName = textFileName + "_" + string(trellisPolynomial(polyIter));
            end
            textFileName = textFileName + "_CRC_";
            textFileName = textFileName + string(crcPolynomial);
            textFileName = textFileName + "_SNR_";
            textFileName = textFileName + sprintf('%.2f',round(snrDBList(snrIter),2));
            textFileName = textFileName + "_TRIALS_";
            textFileName = textFileName + string(numOfTrials);
            textFileName = textFileName + "_K_";
            textFileName = textFileName + string(k);
            textFileName = textFileName + "_LIST_";
            textFileName = textFileName + string(maxListSize(iter0));
            textFileName = textFileName + "_ERRORS_";
            textFileName = textFileName + string(minNumberOfErrors);
            textFileName = textFileName + "_TBCC_";
            textFileName = textFileName + string(tailBitingFlag);
            textFileName = textFileName + "_PUNC_";
            textFileName = textFileName + string(punctureFlag);

            if (punctureFlag == 1)
                textFileName = textFileName + "_" + string(punctureRepetitionLength);
                textFileName = textFileName + "_" + string(punctureCount);
                for iter2 = 1:length(punctureNodeLocations)
                    textFileName = textFileName + "_" + string(punctureNodeLocations(iter2));
                end
            end
            
            if (~(inputSequenceLength == 0))
                textFileName = textFileName + "_SEQUENCE_";
                for iterSequence = 1:inputSequenceLength
                    textFileName = textFileName + string(inputSequence(iterSequence));
                end
            end

            textFileName = textFileName + ".txt";

            fileID = fopen(textFileName, 'r');
            firstLine = fscanf(fileID, "%d %d %d %f %d %d %d %d %d %d %f %f", [1 12]);
            %theRest = fscanf(fileID, "%d %d %d %d\n",[maxListSize(iter0) 4]);
            theRest = fscanf(fileID, "%d");
            sim_data.num_of_trials = firstLine(1);
            sim_data.k = firstLine(2);
            sim_data.n = firstLine(3);
            sim_data.snr = firstLine(4);
            sim_data.num_bit_errors = firstLine(5);
            sim_data.num_frame_errors = firstLine(6);
            sim_data.num_crc_check = firstLine(7);
            sim_data.largest_list_depth = firstLine(8);
            sim_data.max_list_size = firstLine(9);
            sim_data.total_insertions = firstLine(10);
            sim_data.time_trellis = firstLine(11);
            sim_data.time_traceback = firstLine(12);

            sim_data.listDepth = theRest(1:4:end,:);
            sim_data.num_ue = theRest(2:4:end,:);
            sim_data.num_nack = theRest(3:4:end,:);
            sim_data.num_fe = theRest(4:4:end,:);
            fclose(fileID);

            % Saving the data to a final location
            dataFileName = "T";
            for polyIter = 1:length(trellisPolynomial)
                dataFileName = dataFileName + "_" + string(trellisPolynomial(polyIter));
            end
            dataFileName = dataFileName + "_CRC_";
            dataFileName = dataFileName + string(crcPolynomial);
            dataFileName = dataFileName + "_SNR_";
            dataFileName = dataFileName + sprintf('%.2f',round(snrDBList(snrIter),2));
            %dataFileName = dataFileName + "_";
            %dataFileName = dataFileName + sprintf('%.2f',round(snrDBList(snrIter),2));
            dataFileName = dataFileName + "_TRIALS_";
            dataFileName = dataFileName + string(numOfTrials);
            dataFileName = dataFileName + "_K_";
            dataFileName = dataFileName + string(k);
            dataFileName = dataFileName + "_LIST_";
            dataFileName = dataFileName + string(maxListSize(iter0));
            dataFileName = dataFileName + "_";
            dataFileName = dataFileName + string(maxListSize(length(maxListSize)));
            dataFileName = dataFileName + "_ERRORS_";
            dataFileName = dataFileName + string(minNumberOfErrors);

            if (~(inputSequenceLength == 0))
                dataFileName = dataFileName + "_SEQUENCE_";
                for iterSequence = 1:inputSequenceLengths
                    dataFileName = dataFileName + string(inputSequence(iterSequence));
                end
            end

            dataFileName = dataFileName + ".mat";

            save(dataFileName);
        end
    end
end